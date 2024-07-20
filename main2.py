import ezdxf


def read_dwg_and_filter_polylines(file_path):
    # Read the DWG file
    doc = ezdxf.readfile(file_path)

    # Get the model space
    msp = doc.modelspace()

    # Filter closed polylines
    closed_polylines = [entity for entity in msp.query('LWPOLYLINE') if entity.closed]

    return closed_polylines


def calculate_boundaries(polylines):
    min_x = min_y = float('inf')
    max_x = max_y = float('-inf')

    for polyline in polylines:
        for point in polyline.get_points():
            x, y = point[:2]  # Get only x and y coordinates
            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)

    # Add some padding (e.g., 10% of the range)
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


def write_section_info(polylines, output_file):
    boundaries = calculate_boundaries(polylines)

    with open(output_file, 'w') as f:
        f.write("# Begin section definition.\n")
        f.write("Begin_Section\n")
        f.write("\tBegin_Builder\n")
        f.write("\t\tNAME = Section1\n\n")

        # Write builder properties
        f.write(f"\t\tBoundary_Bottom = {boundaries['Boundary_Bottom']:.2f}\n")
        f.write(f"\t\tBoundary_Left = {boundaries['Boundary_Left']:.2f}\n")
        f.write(f"\t\tBoundary_Right = {boundaries['Boundary_Right']:.2f}\n")
        f.write(f"\t\tBoundary_Top = {boundaries['Boundary_Top']:.2f}\n")
        f.write("\t\tMin_Triangle_Area = 64.52\n")  # This might need adjustment
        f.write("\t\tMax_Number_of_Fibers = 4000\n\n")  # This might need adjustment

        # Write window properties
        f.write(f"\t\tWindow_Left = {boundaries['Window_Left']:.2f}\n")
        f.write(f"\t\tWindow_Bottom = {boundaries['Window_Bottom']:.2f}\n")
        f.write(f"\t\tWindow_Height = {boundaries['Window_Height']:.2f}\n")
        f.write("\tEnd_Builder\n")

        f.write("\tBegin_UserComments\n")
        f.write("\tEnd_UserComments\n")

        # Write shapes
        for i, polyline in enumerate(polylines):
            f.write("\tBegin_Shape\n")
            f.write(f"\t\tMATERIAL = BiLinear{i + 1}\n")
            f.write("\t\tMESH = 47.75\n")
            f.write("\t\tBegin_Line\n")

            # Write polyline points
            for point in polyline.get_points():
                f.write(f"\t\t\t{point[0]:.1f}, {point[1]:.1f}\n")

            # Write the first point again to close the shape
            first_point = polyline.get_points()[0]
            f.write(f"\t\t\t{first_point[0]:.1f}, {first_point[1]:.1f}\n")

            f.write("\t\tEnd_Line\n")
            f.write("\tEnd_Shape\n")

        f.write("End_Section\n")


# Usage
dwg_file = "section5.dxf"
output_file = "section_info.txt"

closed_polylines = read_dwg_and_filter_polylines(dwg_file)
write_section_info(closed_polylines, output_file)