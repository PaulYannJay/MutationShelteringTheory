import pandas as pd
import argparse

def main(input_file, output_file):
    # Define the file path and the column names
    file_path = input_file
    col_names = ["N", "mu", "r", "h", "s", "sInv", "sim_cycle", "start", "end", "Rep", "Chrom", 
                 "meanNbMutWInv", "minNbMutWInv", "maxNbMutWInv", "sdNbMutWInv", "meanFreqMutWInv", 
                 "meanNbMutNoInv", "minNbMutNoInv", "maxNbMutNoInv", "sdNbMutNoInv", "meanFreqMutNoInv", 
                 "meanFitnessInv", "meanFitnessNoInv", "Freq"]
    
    # Load the entire dataset at once
    df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=col_names)
    
    # Calculate the "Size" column and ensure 'Rep' is an integer
    df["Size"] = df["end"] - df["start"]
    df["Rep"] = df["Rep"].astype(int)

    # Group and summarize without chunking
    grouped = df.groupby(["N", "mu", "r", "h", "s", "sInv", "Size", "Chrom", "Rep"]).agg(
        MaxFreq=("Freq", "max"),
        MaxGen=("sim_cycle", "max")
    ).reset_index()
    grouped["n"] = 1
    
    # Separate by chromosome type
    dfSumSumAuto = grouped[grouped["Chrom"] == "Autosome"].groupby(["N", "mu", "r", "h", "s", "sInv", "Size", "Chrom"]).agg(
        Fixed=("MaxFreq", lambda x: (x >= 0.95).sum()),
        NotFixed=("MaxFreq", lambda x: (x < 0.95).sum()),
        n=("n", "sum")
    ).reset_index()
    
    dfSumSumSex = grouped[grouped["Chrom"] == "SexChrom"].groupby(["N", "mu", "r", "h", "s", "sInv", "Size", "Chrom"]).agg(
        Fixed=("MaxFreq", lambda x: (x >= 0.25).sum()),
        NotFixed=("MaxFreq", lambda x: (x < 0.25).sum()),
        n=("n", "sum")
    ).reset_index()
    
    # Modify the Chromosome labels
    dfSumSumAuto["Chrom"] = "Autosome"
    dfSumSumSex["Chrom"] = "Y chromosome"
    
    # Combine both dataframes
    dfSumSumAll = pd.concat([dfSumSumAuto, dfSumSumSex])
    
    # Output the final aggregated table
    dfSumSumAll.reset_index(drop=True, inplace=True)
    print(dfSumSumAll)
    
    # Optionally, save the final result to a CSV file
    dfSumSumAll.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summmarise")
    parser.add_argument('-o', '--output', required=True, help="File path for output")
    parser.add_argument('-i', '--input', required=True, help="File path for input")
    
    args = parser.parse_args()
    
    main(args.input, args.output)

