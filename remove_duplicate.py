import sys
import pandas as pd

def main():
    counts_f = sys.argv[1]
    out_f = sys.argv[2]

    df = pd.read_csv(counts_f, sep='\t', index_col=0)
    print(df)
    df = df.drop(labels=['C054'], axis=1)
    print(df)
    df.to_csv(out_f, sep='\t')

if __name__ == "__main__":
    main()
