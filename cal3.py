import os
import matplotlib.pyplot as plt
import ezdxf
import numpy as np
from shapely.geometry import Polygon
from shapely.geometry import Point
from tqdm import tqdm
import triangle


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection



from matplotlib.patches import Polygon as MplPolygon
from matplotlib.collections import PatchCollection


def create_mesh(outer_polyline, inner_polylines, max_area):
    # 准备输入数据
    vertices = get_coordinates(outer_polyline)
    segments = []
    holes = []

    # 添加外部轮廓
    for i in range(len(vertices)):
        segments.append([i, (i + 1) % len(vertices)])

    # 添加内部轮廓和孔
    offset = len(vertices)
    for inner in inner_polylines:
        inner_vertices = get_coordinates(inner)
        holes.append(Polygon(inner_vertices).representative_point().coords[0])
        for i in range(len(inner_vertices)):
            vertices.append(inner_vertices[i])
            segments.append([offset + i, offset + (i + 1) % len(inner_vertices)])
        offset += len(inner_vertices)

    # 创建三角剖分
    mesh = triangle.triangulate({'vertices': vertices, 'segments': segments, 'holes': holes},
                                f'pqa{max_area}')

    print(f"网格三角形数量: {len(mesh['triangles'])}")
    print(f"网格顶点数量: {len(mesh['vertices'])}")

    return mesh


def get_coordinates(polyline):
    return [(vertex[0], vertex[1]) for vertex in polyline]
def create_grid(polyline, grid_size):
    min_x, max_x, min_y, max_y = calculate_boundaries(polyline)
    x = np.arange(min_x, max_x, grid_size)
    y = np.arange(min_y, max_y, grid_size)
    return np.meshgrid(x, y)



def is_point_inside(point, outer_polyline, inner_polylines):
    point = Point(point)
    outer_poly = Polygon(get_coordinates(outer_polyline))
    if not outer_poly.contains(point):
        return False
    for inner in inner_polylines:
        inner_poly = Polygon(get_coordinates(inner))
        if inner_poly.contains(point):
            return False
    return True


def calculate_moments_of_inertia_mesh(mesh, outer_polyline, inner_polylines, centroid):
    cx, cy = centroid
    Ix = Iy = Ixy = 0
    area = 0

    outer_polygon = Polygon(get_coordinates(outer_polyline))
    inner_polygons = [Polygon(get_coordinates(inner)) for inner in inner_polylines]

    total_triangles = len(mesh['triangles'])
    progress_bar = tqdm(total=total_triangles, desc="计算惯性矩", unit="三角形")

    for triangle in mesh['triangles']:
        p1, p2, p3 = [mesh['vertices'][i] for i in triangle]
        triangle_centroid = [(p1[0] + p2[0] + p3[0]) / 3, (p1[1] + p2[1] + p3[1]) / 3]

        # 检查三角形是否在截面内（在外部多边形内但不在任何内部多边形内）
        if outer_polygon.contains(Point(triangle_centroid)) and not any(
                inner.contains(Point(triangle_centroid)) for inner in inner_polygons):
            a = triangle_area(p1, p2, p3)
            x, y = triangle_centroid

            dx, dy = x - cx, y - cy
            Ix += (dy ** 2 * a)
            Iy += (dx ** 2 * a)
            Ixy += (dx * dy * a)
            area += a

        progress_bar.update(1)

    progress_bar.close()

    print(f"网格计算的总面积: {area:.2f}")
    if total_triangles > 0:
        print(
            f"最小三角形面积: {min(triangle_area(mesh['vertices'][i], mesh['vertices'][j], mesh['vertices'][k]) for i, j, k in mesh['triangles']):.6f}")
        print(
            f"最大三角形面积: {max(triangle_area(mesh['vertices'][i], mesh['vertices'][j], mesh['vertices'][k]) for i, j, k in mesh['triangles']):.6f}")

    return Ix, Iy, Ixy, area

def triangle_area(p1, p2, p3):
    return 0.5 * abs((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p3[0] - p1[0]) * (p2[1] - p1[1]))
def read_dwg_and_filter_polylines(file_path):
    doc = ezdxf.readfile(file_path)
    msp = doc.modelspace()
    closed_polylines = [entity for entity in msp.query('LWPOLYLINE') if entity.closed]
    return closed_polylines


def calculate_boundaries(polyline):
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')
    for point in polyline.get_points():
        x, y = point[:2]
        min_x = min(min_x, x)
        max_x = max(max_x, x)
        min_y = min(min_y, y)
        max_y = max(max_y, y)
    return min_x, max_x, min_y, max_y


def get_outermost_polyline(polylines):
    def bounding_area(polyline):
        min_x, max_x, min_y, max_y = calculate_boundaries(polyline)
        return (max_x - min_x) * (max_y - min_y)

    return max(polylines, key=bounding_area)


def is_point_inside_polyline(point, polyline):
    x, y = point
    points = polyline.get_points()
    n = len(points)
    inside = False
    p1x, p1y = points[0][:2]
    for i in range(n + 1):
        p2x, p2y = points[i % n][:2]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def classify_polylines(polylines):
    outermost_polyline = get_outermost_polyline(polylines)
    inner_polylines = [polyline for polyline in polylines if
                       polyline is not outermost_polyline and is_point_inside_polyline(polyline.get_points()[0][:2],
                                                                                       outermost_polyline)]
    return outermost_polyline, inner_polylines

def calculate_area_numerical(outer_polyline, inner_polylines, grid_size):
    x_grid, y_grid = create_grid(outer_polyline, grid_size)
    area = 0
    for i in range(x_grid.shape[0]):
        for j in range(x_grid.shape[1]):
            if is_point_inside((x_grid[i, j], y_grid[i, j]), outer_polyline, inner_polylines):
                area += grid_size**2
    return area

def calculate_area_and_centroid(polyline):
    points = polyline.get_points()
    n = len(points)
    area = 0
    cx = cy = 0
    for i in range(n):
        x0, y0 = points[i][:2]
        x1, y1 = points[(i + 1) % n][:2]
        a = x0 * y1 - x1 * y0
        area += a
        cx += (x0 + x1) * a
        cy += (y0 + y1) * a
    area *= 0.5
    cx = cx / (6 * area)
    cy = cy / (6 * area)
    return abs(area), (cx, cy)


def calculate_moments_of_inertia(outer_polyline, inner_polylines, centroid):
    def polygon_moments_of_inertia(polygon, centroid):
        cx, cy = centroid
        Ix = Iy = Ixy = 0
        for i in range(len(polygon.exterior.coords) - 1):
            x0, y0 = polygon.exterior.coords[i]
            x1, y1 = polygon.exterior.coords[i + 1]
            x0 -= cx
            y0 -= cy
            x1 -= cx
            y1 -= cy
            common = x0 * y1 - x1 * y0
            Ix += (y0 ** 2 + y0 * y1 + y1 ** 2) * common
            Iy += (x0 ** 2 + x0 * x1 + x1 ** 2) * common
            Ixy += (x0 * y1 + 2 * x0 * y0 + 2 * x1 * y1 + x1 * y0) * common
        Ix /= 12
        Iy /= 12
        Ixy /= 24
        return abs(Ix), abs(Iy), abs(Ixy)

    def get_coordinates(polyline):
        return [(x, y) for x, y, *_ in polyline.get_points()]

    outer_polygon = Polygon(get_coordinates(outer_polyline))
    Ix_outer, Iy_outer, Ixy_outer = polygon_moments_of_inertia(outer_polygon, centroid)

    Ix_inner_total = Iy_inner_total = Ixy_inner_total = 0
    for inner_polyline in inner_polylines:
        inner_polygon = Polygon(get_coordinates(inner_polyline))
        Ix_inner, Iy_inner, Ixy_inner = polygon_moments_of_inertia(inner_polygon, centroid)
        Ix_inner_total += Ix_inner
        Iy_inner_total += Iy_inner
        Ixy_inner_total += Ixy_inner

    Ix_total = Ix_outer - Ix_inner_total
    Iy_total = Iy_outer - Iy_inner_total
    Ixy_total = Ixy_outer - Ixy_inner_total
    return Ix_total, Iy_total, Ixy_total


def calculate_polar_moment_of_inertia(Ix, Iy, Ixy):
    return Ix + Iy
def calculate_total_area(outer_polyline, inner_polylines):
    outer_area, _ = calculate_area_and_centroid(outer_polyline)
    inner_areas = [calculate_area_and_centroid(inner)[0] for inner in inner_polylines]
    return outer_area - sum(inner_areas)

def draw_section(outer_polyline, inner_polylines):
    def plot_polyline(ax, polyline, color='blue', fill_color=None):
        points = polyline.get_points()
        x, y = zip(*[(point[0], point[1]) for point in points])
        ax.plot(x, y, color=color)
        if fill_color:
            ax.fill(x, y, color=fill_color, alpha=0.3)

    fig, ax = plt.subplots()
    plot_polyline(ax, outer_polyline, 'blue', 'blue')
    for inner_polyline in inner_polylines:
        plot_polyline(ax, inner_polyline, 'red', 'white')

    ax.set_aspect('equal')
    plt.show()


def visualize_grid_points(outer_polyline, inner_polylines, grid_size):
    x_grid, y_grid = create_grid(outer_polyline, grid_size)
    inside_points = []
    for i in range(x_grid.shape[0]):
        for j in range(x_grid.shape[1]):
            if is_point_inside((x_grid[i, j], y_grid[i, j]), outer_polyline, inner_polylines):
                inside_points.append((x_grid[i, j], y_grid[i, j]))

    fig, ax = plt.subplots()
    draw_section(outer_polyline, inner_polylines)
    x, y = zip(*inside_points)
    ax.scatter(x, y, color='red', s=1)
    plt.show()




def visualize_mesh(mesh, outer_polyline, inner_polylines, output_file):
    fig, ax = plt.subplots(figsize=(12, 12))

    # 调试信息
    print("Outer polyline shape:", np.array(outer_polyline).shape)
    print("First few points of outer polyline:", outer_polyline[:5])

    # 绘制外部轮廓
    outer_polyline = np.array(outer_polyline)
    if outer_polyline.ndim == 2 and outer_polyline.shape[1] == 2:
        ax.plot(outer_polyline[:, 0], outer_polyline[:, 1], 'k-', linewidth=2)
    elif outer_polyline.ndim == 1 and len(outer_polyline) % 2 == 0:
        outer_polyline = outer_polyline.reshape(-1, 2)
        ax.plot(outer_polyline[:, 0], outer_polyline[:, 1], 'k-', linewidth=2)
    else:
        print("Error: Unexpected format for outer_polyline")
        return

    # 绘制内部轮廓
    for inner in inner_polylines:
        inner = np.array(inner)
        if inner.ndim == 2 and inner.shape[1] == 2:
            ax.plot(inner[:, 0], inner[:, 1], 'k-', linewidth=2)
        elif inner.ndim == 1 and len(inner) % 2 == 0:
            inner = inner.reshape(-1, 2)
            ax.plot(inner[:, 0], inner[:, 1], 'k-', linewidth=2)
        else:
            print(f"Error: Unexpected format for inner polyline: {inner.shape}")

    # 绘制网格边缘
    edges = set()
    for triangle in mesh['triangles']:
        for i in range(3):
            edge = tuple(sorted([triangle[i], triangle[(i + 1) % 3]]))
            edges.add(edge)

    lines = [[mesh['vertices'][i], mesh['vertices'][j]] for i, j in edges]
    lc = LineCollection(lines, colors='gray', linewidths=0.5, alpha=0.5)
    ax.add_collection(lc)

    # 设置轴的范围和外观
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Mesh Visualization')

    # 保存图像
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

    print(f"Mesh visualization saved to {output_file}")
def debug_geometry(outer_polyline, inner_polylines):
    print("外部多边形坐标:")
    for point in outer_polyline:
        print(f"  ({point[0]:.2f}, {point[1]:.2f})")

    print("\n内部多边形坐标:")
    for i, inner in enumerate(inner_polylines):
        print(f"内部多边形 {i + 1}:")
        for point in inner:
            print(f"  ({point[0]:.2f}, {point[1]:.2f})")

    min_x, max_x, min_y, max_y = calculate_boundaries(outer_polyline)
    print(f"\n坐标范围: X({min_x:.2f}, {max_x:.2f}), Y({min_y:.2f}, {max_y:.2f})")
    print(f"外部多边形周长: {calculate_perimeter(outer_polyline):.2f}")


def calculate_perimeter(polyline):
    return sum(
        ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5 for p1, p2 in zip(polyline, polyline[1:] + [polyline[0]]))

if __name__ == "__main__":
    dwg_file = "section5.dxf"
    grid_size = 5  # 根据需要调整网格大小

    closed_polylines = read_dwg_and_filter_polylines(dwg_file)
    outer_polyline, inner_polylines = classify_polylines(closed_polylines)

    print("计算初始重心...")
    outer_area, outer_centroid = calculate_area_and_centroid(outer_polyline)
    inner_areas_and_centroids = [calculate_area_and_centroid(inner_polyline) for inner_polyline in inner_polylines]
    total_inner_area = sum(area for area, _ in inner_areas_and_centroids)
    total_area = outer_area - total_inner_area
    cx_total = (outer_area * outer_centroid[0] - sum(
        area * centroid[0] for area, centroid in inner_areas_and_centroids)) / total_area
    cy_total = (outer_area * outer_centroid[1] - sum(
        area * centroid[1] for area, centroid in inner_areas_and_centroids)) / total_area
    centroid_total = (cx_total, cy_total)

    print("计算围合面积...")
    enclosed_area = calculate_total_area(outer_polyline, inner_polylines)

    print("开始数值积分计算...")
    print("生成网格...")
    max_area = 100  # 控制网格的细密程度，可以根据需要调整
    mesh = create_mesh(outer_polyline, inner_polylines, max_area)

    print("计算惯性矩...")
    Ix, Iy, Ixy, area_mesh = calculate_moments_of_inertia_mesh(mesh, outer_polyline, inner_polylines, centroid_total)
    J = Ix + Iy

    print("\n计算结果：")
    print(f"围合面积: {enclosed_area:.2f}")
    print(f"网格化面积: {area_mesh:.2f}")
    print(f"面积差异: {abs(enclosed_area - area_mesh):.2f} ({abs(enclosed_area - area_mesh) / enclosed_area * 100:.2f}%)")
    print(f"重心: ({centroid_total[0]:.2f}, {centroid_total[1]:.2f})")
    print(f"绕 x 轴的抗弯惯性矩 (Ix): {Ix:.2f}")
    print(f"绕 y 轴的抗弯惯性矩 (Iy): {Iy:.2f}")
    print(f"惯性积 (Ixy): {Ixy:.2f}")
    print(f"极惯性矩 (J): {J:.2f}")

    # 添加额外的检查
    print("\n额外检查：")
    print(f"Ix/area: {Ix/area_mesh:.2f}")
    print(f"Iy/area: {Iy/area_mesh:.2f}")
    print(f"J/area: {J/area_mesh:.2f}")

    # 可视化网格
    visualize_mesh(mesh, outer_polyline, inner_polylines, "mesh_visualization.png")