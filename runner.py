from pipeline import run_pgx_technical_report

def main():
    sample = "EDX2508083837"
    gene_phenotypes = run_pgx_technical_report(sample, return_all_steps=True)
    print(" Gene Phenotypes : ",gene_phenotypes)

if __name__ == "__main__":
    main()