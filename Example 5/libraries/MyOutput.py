import csv


def saveCSV(table, filename):
    filename = "output/" + filename
    s = str(table)

    result = [
        tuple(filter(None, map(str.strip, splitline)))
        for line in s.splitlines()
        for splitline in [line.split("|")]
        if len(splitline) > 1
    ]

    with open(filename, "w", encoding="utf-8") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerows(result)