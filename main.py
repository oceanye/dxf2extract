import ezdxf


def read_dwg_and_filter_polylines(file_path):
    # Read the DWG file
    doc = ezdxf.readfile(file_path)

    # Get the model space
    msp = doc.modelspace()

    # Filter closed polylines
    closed_polylines = [entity for entity in msp.query('LWPOLYLINE') if entity.closed]

    return closed_polylines


def write_section_info(polylines, output_file):
    with open(output_file, 'w') as f:
        f.write("# Begin section definition.\n")
        f.write("Begin_Section\n")
        f.write("\tBegin_Builder\n")
        f.write("\t\tNAME = Section1\n\n")

        # Write builder properties (you may need to adjust these values)
        f.write("\t\tBoundary_Bottom = -2540\n")
        f.write("\t\tBoundary_Left = -2540\n")
        f.write("\t\tBoundary_Right = 2540\n")
        f.write("\t\tBoundary_Top = 2540\n")
        f.write("\t\tMin_Triangle_Area = 64.52\n")
        f.write("\t\tMax_Number_of_Fibers = 4000\n\n")

        # Write window properties (you may need to adjust these values)
        f.write("\t\tWindow_Left = 304.8\n")
        f.write("\t\tWindow_Bottom = 304.8\n")
        f.write("\t\tWindow_Height = 1219\n")
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
dwg_file = "section4.dxf"
output_file = "section_info.txt"

closed_polylines = read_dwg_and_filter_polylines(dwg_file)
write_section_info(closed_polylines, output_file)