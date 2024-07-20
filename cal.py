import ezdxf
from shapely.geometry import Polygon, MultiPolygon, Point
import numpy as np

def read_regions_from_dxf(file_path):
    """从DXF文件中读取所有region"""
    doc = ezdxf.readfile(file_path)
    msp = doc.modelspace()

    regions = []
    for entity in msp:
        if entity.dxftype() == 'HATCH':
            for path in entity.paths:
                if path.is_closed:
                    polygon = Polygon(path.vertices)
                    regions.append(polygon)
                else:
                    print("非闭合路径被忽略")

    return MultiPolygon(regions)

def calculate_moment_of_inertia(shape):
    """计算惯性矩"""
    if shape.geom_type == 'MultiPolygon':
        Ix, Iy = 0, 0
        for polygon in shape.geoms:
            ix, iy = calculate_polygon_moment_of_inertia(polygon)
            Ix += ix
            Iy += iy
    else:
        Ix, Iy = calculate_polygon_moment_of_inertia(shape)

    return Ix, Iy

def calculate_polygon_moment_of_inertia(polygon):
    """计算单个多边形的惯性矩"""
    coords = np.array(polygon.exterior.coords)
    x, y = coords[:, 0], coords[:, 1]

    # 使用Green定理计算
    Ix = (1 / 12) * np.sum((y[:-1] ** 2 + y[:-1] * y[1:] + y[1:] ** 2) * (x[:-1] * y[1:] - x[1:] * y[:-1]))
    Iy = (1 / 12) * np.sum((x[:-1] ** 2 + x[:-1] * x[1:] + x[1:] ** 2) * (x[:-1] * y[1:] - x[1:] * y[:-1]))

    # 减去内部孔洞的惯性矩
    for interior in polygon.interiors:
        coords = np.array(interior.coords)
        x, y = coords[:, 0], coords[:, 1]
        Ix -= (1 / 12) * np.sum((y[:-1] ** 2 + y[:-1] * y[1:] + y[1:] ** 2) * (x[:-1] * y[1:] - x[1:] * y[:-1]))
        Iy -= (1 / 12) * np.sum((x[:-1] ** 2 + x[:-1] * x[1:] + x[1:] ** 2) * (x[:-1] * y[1:] - x[1:] * y[:-1]))

    return abs(Ix), abs(Iy)

def calculate_area_and_centroid(shape):
    """计算区域的面积和重心"""
    if shape.geom_type == 'MultiPolygon':
        total_area = 0
        centroid_x, centroid_y = 0, 0
        for polygon in shape.geoms:
            area = polygon.area
            centroid = polygon.centroid
            total_area += area
            centroid_x += area * centroid.x
            centroid_y += area * centroid.y
        centroid_x /= total_area
        centroid_y /= total_area
    else:
        total_area = shape.area
        centroid = shape.centroid
        centroid_x, centroid_y = centroid.x, centroid.y

    return total_area, Point(centroid_x, centroid_y)

def main(dxf_file_path):
    # 读取DXF文件中的region
    shape = read_regions_from_dxf(dxf_file_path)

    # 计算惯性矩
    Ix, Iy = calculate_moment_of_inertia(shape)
    print(f"惯性矩 Ix = {Ix}")
    print(f"惯性矩 Iy = {Iy}")

    # 计算面积和重心
    area, centroid = calculate_area_and_centroid(shape)
    print(f"面积 = {area}")
    print(f"重心 = ({centroid.x}, {centroid.y})")

if __name__ == "__main__":
    dxf_file_path = "section.dxf"  # 替换为你的DXF文件路径
    main(dxf_file_path)