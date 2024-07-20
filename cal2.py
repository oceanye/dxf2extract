import os
import matplotlib.pyplot as plt
import ezdxf
import numpy as np
from shapely.geometry import Polygon


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


# 使用示例
dwg_file = "section5.dxf"

closed_polylines = read_dwg_and_filter_polylines(dwg_file)
outer_polyline, inner_polylines = classify_polylines(closed_polylines)
outer_area, outer_centroid = calculate_area_and_centroid(outer_polyline)
inner_areas_and_centroids = [calculate_area_and_centroid(inner_polyline) for inner_polyline in inner_polylines]
total_inner_area = sum(area for area, _ in inner_areas_and_centroids)

# 计算总面积和重心
total_area = outer_area - total_inner_area
cx_total = (outer_area * outer_centroid[0] - sum(
    area * centroid[0] for area, centroid in inner_areas_and_centroids)) / total_area
cy_total = (outer_area * outer_centroid[1] - sum(
    area * centroid[1] for area, centroid in inner_areas_and_centroids)) / total_area
centroid_total = (cx_total, cy_total)

# 计算惯性矩
Ix, Iy, Ixy = calculate_moments_of_inertia(outer_polyline, inner_polylines, centroid_total)
J = calculate_polar_moment_of_inertia(Ix, Iy, Ixy)

# 输出基本参数
print(f"总面积: {total_area:.2f}")
print(f"重心: ({centroid_total[0]:.2f}, {centroid_total[1]:.2f})")
print(f"绕 x 轴的抗弯惯性矩 (Ix): {Ix:.2f}")
print(f"绕 y 轴的抗弯惯性矩 (Iy): {Iy:.2f}")
print(f"抗扭惯性矩 (Polar Moment of Inertia, J): {J:.2f}")

# 绘制截面
draw_section(outer_polyline, inner_polylines)