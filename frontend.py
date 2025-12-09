"""
Streamlit Frontend for PGX Technical Report
Displays all sub-function outputs and final gene phenotypes in table format
"""

import streamlit as st
import pandas as pd
from pipeline import run_pgx_technical_report
import json

# Page configuration
st.set_page_config(
    page_title="PGX Technical Report",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #ffffff;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.5rem;
        font-weight: bold;
        color: #2c3e50;
        margin-top: 2rem;
        margin-bottom: 1rem;
        padding-bottom: 0.5rem;
        border-bottom: 2px solid #3498db;
    }
    .gene-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 8px;
        margin-bottom: 1rem;
        border-left: 4px solid #3498db;
        color: #000000;
    }
    .gene-card h4 {
        color: #000000;
        margin: 0;
    }
    h3 {
        color: #000000;
    }
    [data-testid="stSidebar"] h3 {
        color: #ffffff !important;
    }
    .success-box {
        background-color: #d4edda;
        padding: 1rem;
        border-radius: 8px;
        border-left: 4px solid #28a745;
        margin-bottom: 1rem;
        color: #000000;
    }
    </style>
""", unsafe_allow_html=True)

def display_filtered_data(filtered_data):
    """Display the filtered PGX data"""
    if filtered_data is not None and not filtered_data.empty:
        st.markdown('<div class="section-header">📊 Filtered PGX Data</div>', unsafe_allow_html=True)
        st.dataframe(filtered_data, use_container_width=True)
        st.info(f"Total rows: {len(filtered_data)}")
    else:
        st.warning("No filtered data available")

def display_unique_genes(unique_gene_list):
    """Display unique gene list"""
    if unique_gene_list:
        st.markdown('<div class="section-header">🧬 Unique Genes</div>', unsafe_allow_html=True)
        genes_df = pd.DataFrame({"Gene": unique_gene_list})
        st.dataframe(genes_df, use_container_width=True)
        st.info(f"Total genes: {len(unique_gene_list)}")
    else:
        st.warning("No genes found")

def display_gene_summary(gene_summary):
    """Display gene summary with SNP calls"""
    if gene_summary:
        st.markdown('<div class="section-header">📋 Gene Summary</div>', unsafe_allow_html=True)
        for gene, info in gene_summary.items():
            with st.expander(f"Gene: {gene}", expanded=False):
                calls = info.get("calls", [])
                if calls:
                    calls_df = pd.DataFrame(calls)
                    st.dataframe(calls_df, use_container_width=True)
                    st.info(f"Total SNP calls: {len(calls)}")
                else:
                    st.warning("No calls found for this gene")
                
def display_star_alleles(star_alleles):
    """Display star alleles mapping"""
    if star_alleles:
        st.markdown('<div class="section-header">⭐ Star Alleles</div>', unsafe_allow_html=True)
        for gene, entries in star_alleles.items():
            if entries:
                with st.expander(f"Gene: {gene} - Star Allele Mappings", expanded=False):
                    star_df = pd.DataFrame(entries)
                    st.dataframe(star_df, use_container_width=True)
                    st.info(f"Total mappings: {len(entries)}")
            else:
                st.warning(f"No star allele mappings found for {gene}")

def display_star_functions(star_functions):
    """Display star allele functions"""
    if star_functions:
        st.markdown('<div class="section-header">⚙️ Star Allele Functions</div>', unsafe_allow_html=True)
        for gene, functions in star_functions.items():
            if functions:
                with st.expander(f"Gene: {gene} - Allele Functions", expanded=False):
                    func_df = pd.DataFrame([
                        {"Star Allele": star, "Function": func}
                        for star, func in functions.items()
                    ])
                    st.dataframe(func_df, use_container_width=True)
            else:
                st.warning(f"No function data found for {gene}")

def display_gene_phenotypes(gene_phenotypes):
    """Display final gene phenotypes in table format"""
    if gene_phenotypes:
        st.markdown('<div class="section-header">🎯 Final Gene Phenotypes</div>', unsafe_allow_html=True)
        
        # Create a summary table for all genes
        summary_data = []
        for gene, diplotypes in gene_phenotypes.items():
            for diplotype, info in diplotypes.items():
                summary_data.append({
                    "Gene": gene,
                    "Diplotype": diplotype,
                    "Alleles": ", ".join(info.get("alleles", [])),
                    "Phenotype": info.get("phenotype", "Unknown"),
                    "Activity Score": info.get("activity_score", ""),
                    "EHR Priority": info.get("ehr_priority", ""),
                    "Allele Functions": ", ".join([
                        f"{star}: {func}" 
                        for star, func in info.get("allele_functions", {}).items()
                    ])
                })
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            st.dataframe(summary_df, use_container_width=True)
            
            # Display by gene in separate tables
            st.markdown("### 📊 Detailed View by Gene")
            st.markdown('<style>h3 { color: #000000 !important; }</style>', unsafe_allow_html=True)
            for gene in sorted(gene_phenotypes.keys()):
                st.markdown(f'<div class="gene-card"><h4 style="color: #000000;">🧬 {gene}</h4></div>', unsafe_allow_html=True)
                
                gene_data = []
                for diplotype, info in gene_phenotypes[gene].items():
                    gene_data.append({
                        "Diplotype": diplotype,
                        "Allele 1": info.get("alleles", [""])[0] if len(info.get("alleles", [])) > 0 else "",
                        "Allele 2": info.get("alleles", [""])[1] if len(info.get("alleles", [])) > 1 else "",
                        "Phenotype": info.get("phenotype", "Unknown"),
                        "Activity Score": info.get("activity_score", ""),
                        "EHR Priority": info.get("ehr_priority", ""),
                        "Allele 1 Function": info.get("allele_functions", {}).get(
                            info.get("alleles", [""])[0] if len(info.get("alleles", [])) > 0 else "", 
                            "Unknown"
                        ),
                        "Allele 2 Function": info.get("allele_functions", {}).get(
                            info.get("alleles", [""])[1] if len(info.get("alleles", [])) > 1 else "", 
                            "Unknown"
                        )
                    })
                
                if gene_data:
                    gene_df = pd.DataFrame(gene_data)
                    st.dataframe(gene_df, use_container_width=True)
                    st.markdown("---")
        else:
            st.warning("No phenotype data available")
    else:
        st.warning("No gene phenotypes found")

def main():
    # Initialize session state
    if "results" not in st.session_state:
        st.session_state.results = None
    if "current_sample_id" not in st.session_state:
        st.session_state.current_sample_id = None
    
    # Header
    st.markdown('<div class="main-header">🧬 PGX Technical Report</div>', unsafe_allow_html=True)
    
    # Sidebar for input
    with st.sidebar:
        st.header("⚙️ Configuration")
        sample_id = st.text_input(
            "Sample ID",
            value="EDX2508083837",
            help="Enter the sample ID to process"
        )
        
        st.markdown("---")
        st.markdown('<h3 style="color: #ffffff;">📝 Instructions</h3>', unsafe_allow_html=True)
        st.markdown("""
        1. Enter a Sample ID in the input field
        2. Click 'Run Analysis' to process
        3. View all intermediate results
        4. Check the final gene phenotypes table
        """)
        
        # Clear results button
        if st.session_state.results is not None:
            if st.button("🗑️ Clear Results", use_container_width=True):
                st.session_state.results = None
                st.session_state.current_sample_id = None
                st.rerun()
    
    # Main content area
    run_button = st.button("🚀 Run Analysis", type="primary", use_container_width=True)
    
    if run_button:
        if not sample_id or sample_id.strip() == "":
            st.error("Please enter a valid Sample ID")
        else:
            with st.spinner(f"Processing Sample ID: {sample_id}..."):
                try:
                    # Run the pipeline with all steps
                    results = run_pgx_technical_report(sample_id.strip(), return_all_steps=True)
                    
                    if results:
                        st.session_state.results = results
                        st.session_state.current_sample_id = sample_id.strip()
                        st.rerun()
                    else:
                        st.error("No results returned from analysis")
                        
                except Exception as e:
                    st.error(f"❌ Error during analysis: {str(e)}")
                    st.exception(e)
    
    # Display results if available
    if st.session_state.results is not None:
        results = st.session_state.results
        sample_id = st.session_state.current_sample_id
        
        st.markdown('<div class="success-box">✅ Analysis completed successfully!</div>', unsafe_allow_html=True)
        st.markdown(f"**Sample ID:** {sample_id}")
        
        # Display all sub-function outputs
        st.markdown("## 📊 Analysis Results")
        
        # 1. Filtered Data
        display_filtered_data(results.get("filtered_data"))
        
        # 2. Unique Genes
        display_unique_genes(results.get("unique_gene_list"))
        
        # 3. Gene Summary
        display_gene_summary(results.get("gene_summary"))
        
        # 4. Star Alleles
        display_star_alleles(results.get("star_alleles"))
        
        # 5. Star Functions
        display_star_functions(results.get("star_functions"))
        
        # 6. Final Gene Phenotypes (main output)
        st.markdown("---")
        display_gene_phenotypes(results.get("gene_phenotypes"))
        
        # Download button for results
        st.markdown("---")
        st.markdown("### 💾 Download Results")
        results_json = json.dumps(results.get("gene_phenotypes"), indent=2, default=str)
        st.download_button(
            label="📥 Download Gene Phenotypes (JSON)",
            data=results_json,
            file_name=f"gene_phenotypes_{sample_id}.json",
            mime="application/json"
        )
    
    elif not run_button:
        # Initial state - show instructions
        st.info("👆 Enter a Sample ID and click 'Run Analysis' to begin")
        
        # Show example
        with st.expander("📖 Example Usage", expanded=False):
            st.code("""
Sample ID: EDX2508083837

The analysis will:
1. Filter PGX data by Sample ID
2. Extract unique genes
3. Summarize gene calls
4. Map to star alleles
5. Get star allele functions
6. Map to phenotypes
7. Display final results in tables
            """)

if __name__ == "__main__":
    main()
