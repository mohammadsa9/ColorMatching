import xlsxwriter
import texttable
import os
import shutil
import pathlib


OUTPUT = 0


def makeout():
    global OUTPUT
    OUTPUT = 1
    clear_folder("output")
    clear_folder("output_img")


def clear_folder(name):
    path = os.path.abspath(os.getcwd()) + "\\" + name + "\\"
    try:
        shutil.rmtree(path)
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    except Exception:
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def repeat(input, repeat):
    array = []
    for i in range(repeat):
        array.append(input)
    return array


def save(table, filename):

    look_table = 0
    look_table = texttable.Texttable(100)
    look_table.set_cols_dtype(repeat("t", len(table[0])))
    look_table.set_cols_align(repeat("c", len(table[0])))
    look_table.set_cols_valign(repeat("m", len(table[0])))
    for item in table:
        look_table.add_row(item)
    print(look_table.draw())

    global OUTPUT
    if OUTPUT:
        filename = "output/" + filename + ".xlsx"
        workbook = xlsxwriter.Workbook(filename)
        worksheet = workbook.add_worksheet()
        cell_format = workbook.add_format()
        cell_format.set_align("center")
        cell_format.set_align("vcenter")

        row = 0
        for item in table:
            for col in range(len(item)):
                worksheet.write(row, col, str(item[col]), cell_format)
            row += 1
        workbook.close()

        try:
            import win32com.client as win32

            excel = win32.gencache.EnsureDispatch("Excel.Application")
            wb = excel.Workbooks.Open(os.path.abspath(os.getcwd()) + "\\" + filename)
            ws = wb.Worksheets("Sheet1")
            ws.Columns.AutoFit()
            wb.Save()
            excel.Application.Quit()
        except Exception:
            print("Unable to autofit excel file")