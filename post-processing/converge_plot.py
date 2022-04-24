import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-f', '--file', type=str)
    args = parser.parse_args()
    # df = pd.read_csv("convergePlot.csv", delimiter=",")
    print(args.file)
    df = pd.read_csv(args.file, delimiter=",");
    plt.semilogy(df.iter, df.RSS, '.r', label=r"RSS for $\rho$")
    plt.xlabel(f"iterations"); plt.ylabel(f"RSS")
    plt.title("convergence plot")
    # plt.grid()
    plt.legend()
    plt.savefig(args.file +".png")
    plt.show()

if __name__ == '__main__':
    main()


