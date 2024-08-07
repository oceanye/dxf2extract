import ezdxf


def read_dwg_and_filter_polylines(file_path):
    # 读取DWG文件
    doc = ezdxf.readfile(file_path)

    # 获取模型空间
    msp = doc.modelspace()

    # 过滤闭合的多段线
    closed_polylines = [entity for entity in msp.query('LWPOLYLINE') if entity.closed]

    return closed_polylines


def calculate_boundaries(polylines):
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')

    for polyline in polylines:
        for point in polyline.get_points():
            x, y = point[:2]  # 只获取x和y坐标
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)

    # 添加一些填充（例如，范围的10%）
    padding_x = (max_x - min_x) * 0.1
    padding_y = (max_y - min_y) * 0.1

    return {
        'Boundary_Left': min_x - padding_x,
        'Boundary_Right': max_x + padding_x,
        'Boundary_Bottom': min_y - padding_y,
        'Boundary_Top': max_y + padding_y,
        'Window_Left': min_x,
        'Window_Bottom': min_y,
        'Window_Height': max_y - min_y
    }


def get_outermost_polyline(polylines):
    # 计算每个多段线的包围盒面积
    def bounding_area(polyline):
        min_x = min_y = float('inf')
        max_x = max_y = float('-inf')
        for point in polyline.get_points():
            x, y = point[:2]
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
        return (max_x - min_x) * (max_y - min_y)

    # 找到包围盒面积最大的多段线
    return max(polylines, key=bounding_area)


def write_section_info(polylines, output_file):
    boundaries = calculate_boundaries(polylines)
    outermost_polyline = get_outermost_polyline(polylines)

    with open(output_file, 'w') as f:
        f.write("# Begin section definition.\n")
        f.write("Begin_Section\n")
        f.write("\tBegin_Builder\n")
        f.write("\t\tNAME = Section1\n\n")

        # 写入builder属性
        f.write(f"\t\tBoundary_Bottom = {boundaries['Boundary_Bottom']:.2f}\n")
        f.write(f"\t\tBoundary_Left = {boundaries['Boundary_Left']:.2f}\n")
        f.write(f"\t\tBoundary_Right = {boundaries['Boundary_Right']:.2f}\n")
        f.write(f"\t\tBoundary_Top = {boundaries['Boundary_Top']:.2f}\n")
        f.write("\t\tMin_Triangle_Area = 64.52\n")  # 这个值可能需要调整
        f.write("\t\tMax_Number_of_Fibers = 4000\n\n")  # 这个值可能需要调整

        # 写入窗口属性
        f.write(f"\t\tWindow_Left = {boundaries['Window_Left']:.2f}\n")
        f.write(f"\t\tWindow_Bottom = {boundaries['Window_Bottom']:.2f}\n")
        f.write(f"\t\tWindow_Height = {boundaries['Window_Height']:.2f}\n")
        f.write("\tEnd_Builder\n")

        f.write("\tBegin_UserComments\n")
        f.write("\tEnd_UserComments\n")

        # 写入形状
        for i, polyline in enumerate(polylines):
            f.write("\tBegin_Shape\n")
            if polyline is outermost_polyline:
                f.write("\t\tMATERIAL = Confined1\n")
            else:
                f.write("\t\tMATERIAL = Delete\n")
            f.write("\t\tMESH = 47.75\n")
            f.write("\t\tBegin_Line\n")

            # 写入多段线的点
            for point in polyline.get_points():
                f.write(f"\t\t\t{point[0]:.1f}, {point[1]:.1f}\n")

            # 再次写入第一个点以闭合形状
            first_point = polyline.get_points()[0]
            f.write(f"\t\t\t{first_point[0]:.1f}, {first_point[1]:.1f}\n")

            f.write("\t\tEnd_Line\n")
            f.write("\tEnd_Shape\n")

        f.write("End_Section\n")


# 使用示例
dwg_file = "section6.dxf"
output_file = "section_info6.txt"

closed_polylines = read_dwg_and_filter_polylines(dwg_file)
write_section_info(closed_polylines, output_file)