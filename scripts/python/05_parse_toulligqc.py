#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

FIELD_MAP = {
    "Sample_Name": "Sample",
    "Read_Count": "basecaller.sequencing.summary.1d.extractor.read.count",
    "Read_Pass_Count": "basecaller.sequencing.summary.1d.extractor.read.pass.count",
    "Read_Pass_Percent": "basecaller.sequencing.summary.1d.extractor.read.pass.frequency",
    "Yield": "basecaller.sequencing.summary.1d.extractor.yield",
    "N50": "basecaller.sequencing.summary.1d.extractor.n50",
    "L50": "basecaller.sequencing.summary.1d.extractor.l50",
    "Pass_Read_Mean_Length": "basecaller.sequencing.summary.1d.extractor.pass.reads.sequence.length.mean",
    "Pass_Read_Mean_Length_Min": "basecaller.sequencing.summary.1d.extractor.pass.reads.sequence.length.min",
    "Pass_Read_Mean_Length_Max": "basecaller.sequencing.summary.1d.extractor.pass.reads.sequence.length.max",
    "Pass_Reads_QScore_Mean": "basecaller.sequencing.summary.1d.extractor.pass.reads.mean.qscore.mean",
    "ToulligQC_Version": "toulligqc.info.version"
}

def parse_toulligqc_file(path):
    data = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue
            key, value = line.split('=', 1)
            value = value.strip().strip("[]'")
            try:
                if '.' in value:
                    value = float(value)
                else:
                    value = int(value)
            except ValueError:
                pass
            data[key] = value
    return data

def main():
    parser = argparse.ArgumentParser(description="ToulligQC summary CSV generator")
    parser.add_argument("-i", "--input", required=True, help="Text file with paths to .data/.info files")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    input_file = Path(args.input)
    if not input_file.is_file():
        print(f"Input file {args.input} not found")
        return

    with open(input_file) as f:
        paths = [line.strip() for line in f if line.strip()]

    all_rows = []
    for p in paths:
        file_path = Path(p)
        if not file_path.is_file():
            print(f"Warning: file {p} not found, skipping")
            continue
        info = parse_toulligqc_file(file_path)
        row = {}
        row["Sample_Name"] = file_path.stem
        for csv_col, data_key in FIELD_MAP.items():
            if csv_col == "Sample_Name":
                continue
            row[csv_col] = info.get(data_key, "")
        all_rows.append(row)

    # Write CSV
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(FIELD_MAP.keys()))
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"Wrote {len(all_rows)} samples to {args.output}")

if __name__ == "__main__":
    main()