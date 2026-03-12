import pytest
import pandas as pd
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

from parse import *
from calc import *


class MockArgs:
    """Mock args object for testing"""
    def __init__(self, gct=False, conf=False, redo=False):
        self.gct = gct
        self.conf = conf
        self.redo = redo
        self.cov = None


class TestReadEQTLs:
    """Test suite for reading eQTL files"""
    
    def test_read_eqtls_basic(self):
        """Test basic eQTL file reading"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        
        assert len(eqtl_file.columns) == 19
        assert eqtl_file.columns[0] == "gene_id"
        assert "gene_id_clean" in eqtl_file.columns
        assert "variant_id" in eqtl_file.columns
        
    def test_read_eqtls_gene_id_clean(self):
        """Test that gene versions are properly stripped"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        
        # Check that gene_id_clean doesn't contain version numbers
        assert all('.' not in str(gene_id) for gene_id in eqtl_file.gene_id_clean)
        
    def test_read_eqtls_wrong_format(self):
        """Test that malformed variant IDs raise exception"""
        with pytest.raises(Exception, match="Variant IDs must not begin with a digit"):
            read_eqtls("tests/dataset/test_wrongformat_eqtls.tsv")


class TestVariantIDValidation:
    """Test suite for variant ID format validation"""
    
    def test_is_variant_id_format_valid_correct(self):
        """Test variant ID validation with correct format"""
        eqtl_file = pd.read_csv("tests/dataset/test_eqtls.tsv", sep='\t')
        assert is_variant_id_format_valid(eqtl_file) == False
        
    def test_is_variant_id_format_valid_wrong(self):
        """Test variant ID validation with incorrect format"""
        eqtl_file = pd.read_csv("tests/dataset/test_wrongformat_eqtls.tsv", sep='\t')
        assert is_variant_id_format_valid(eqtl_file) == True


class TestReadExpressions:
    """Test suite for reading expression files"""
    
    def test_read_expressions_basic(self):
        """Test basic expression file reading"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=False)
        
        expr_file = read_expressions("tests/dataset/test_expressions.csv.gz", eqtl_file, args)
        
        assert len(expr_file.columns) == 575  # Name + 574 samples
        assert "Name" in expr_file.columns
        assert not expr_file.empty
        
    def test_read_expressions_gct_format(self):
        """Test reading GCT format expression files"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=True)
        
        expr_file = read_expressions("tests/dataset/test_expressions_gct.csv.gz", eqtl_file, args)
        
        # Description column should be dropped
        assert "Description" not in expr_file.columns
        assert "Name" in expr_file.columns
        
    def test_read_expressions_gene_version_stripping(self):
        """Test that gene versions are stripped from Name column"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=False)
        
        expr_file = read_expressions("tests/dataset/test_expressions.csv.gz", eqtl_file, args)
        
        # Check that Name column doesn't contain version numbers
        assert all('.' not in str(name) for name in expr_file.Name)


class TestDataFiltering:
    """Test suite for data filtering functions"""
    
    def test_drop_noexpr_genes(self):
        """Test dropping genes with no expression data"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=False)
        expr_file = read_expressions("tests/dataset/test_expressions.csv.gz", eqtl_file, args)
        expr_file.index = expr_file.Name
        
        filtered_eqtl = drop_noexpr_genes(eqtl_file, expr_file)
        
        assert filtered_eqtl.columns[0] == 'gene_id'
        assert len(filtered_eqtl) <= len(eqtl_file)
        assert all(gene in expr_file.index for gene in filtered_eqtl.gene_id_clean)


class TestOutputInitialization:
    """Test suite for output dataframe initialization"""
    
    def test_init_output_df_no_conf(self):
        """Test output dataframe initialization without confidence intervals"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(conf=False)
        
        output_df = init_output_df(eqtl_file, args)
        
        assert 'log2_aFC' in output_df.columns
        assert 'log2_aFC_error' in output_df.columns
        assert 'log2_aFC_c0' in output_df.columns
        assert all(pd.isna(output_df.log2_aFC))
        
    def test_init_output_df_with_conf(self):
        """Test output dataframe initialization with confidence intervals"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(conf=True)
        
        output_df = init_output_df(eqtl_file, args)
        
        assert 'log2_aFC' in output_df.columns
        assert 'log2_aFC_min_95_interv' in output_df.columns
        assert 'log2_aFC_plus_95_interv' in output_df.columns
        assert 'log2_aFC_c0_min_95_interv' in output_df.columns
        assert 'log2_aFC_c0_plus_95_interv' in output_df.columns


class TestGeneExpressionNormalization:
    """Test suite for gene expression normalization"""
    
    def test_gene_expression_normalization(self):
        """Test gene expression normalization function"""
        # Create a simple test matrix
        test_data = pd.DataFrame({
            'sample1': [150, 200, 300],
            'sample2': [180, 220, 320],
            'sample3': [160, 210, 310]
        })
        
        normalized = gene_expression_normalization(test_data)
        
        # Check output shape matches input
        assert normalized.shape == test_data.shape
        # Check that normalization happened (values changed)
        assert not np.array_equal(normalized.values, np.log2(test_data + 1).values)


class TestHaplotypeArrays:
    """Test suite for haplotype array manipulation"""
    
    def test_make_haplotype_arrays(self):
        """Test haplotype array creation from VCF format"""
        # Create mock haplotype data in VCF format
        test_hap = pd.DataFrame({
            'sample1': ['0|1:30', '1|0:25'],
            'sample2': ['1|1:35', '0|0:28'],
            'sample3': ['0|0:32', '1|1:30']
        }, index=['variant1', 'variant2'])
        
        h1, h2 = make_haplotype_arrays(test_hap)
        
        # Check shapes
        assert h1.shape == test_hap.shape
        assert h2.shape == test_hap.shape
        
        # Check that values are numeric
        assert h1.dtypes[0] == np.float64
        assert h2.dtypes[0] == np.float64
        
        # Check specific values
        assert h1.loc['variant1', 'sample1'] == 0
        assert h2.loc['variant1', 'sample1'] == 1
        

class TestNaNHandling:
    """Test suite for NaN handling functions"""
    
    def test_get_nan_columns(self):
        """Test identification of columns with NaN values"""
        test_df = pd.DataFrame({
            'col1': [1, 2, 3],
            'col2': [1, np.nan, 3],
            'col3': [1, 2, np.nan],
            'col4': [1, 2, 3]
        })
        
        nan_cols = get_nan_columns(test_df)
        
        assert 'col2' in nan_cols
        assert 'col3' in nan_cols
        assert 'col1' not in nan_cols
        assert 'col4' not in nan_cols
        
    def test_clean_matrix_nans(self):
        """Test cleaning NaN values from matrices"""
        h1 = pd.DataFrame({
            'ind1': [0, 1],
            'ind2': [np.nan, 0],
            'ind3': [1, 1]
        })
        h2 = pd.DataFrame({
            'ind1': [1, 0],
            'ind2': [0, 1],
            'ind3': [0, np.nan]
        })
        expr = pd.DataFrame({
            'ind1': [5.5],
            'ind2': [6.2],
            'ind3': [5.8]
        })
        
        h1_clean, h2_clean, expr_clean = clean_matrix_nans(h1, h2, expr)
        
        # ind2 and ind3 should be dropped
        assert 'ind2' not in h1_clean.columns
        assert 'ind3' not in h2_clean.columns
        assert 'ind1' in h1_clean.columns
        assert len(h1_clean.columns) == 1


class TestLinearEstimate:
    """Test suite for linear estimation functions"""
    
    def test_nocovar_linear_estimate(self):
        """Test linear estimation without covariates"""
        # Create simple test data
        h1 = pd.DataFrame({
            'ind1': [0, 1],
            'ind2': [1, 0],
            'ind3': [0, 0]
        }, index=['var1', 'var2'])
        
        h2 = pd.DataFrame({
            'ind1': [1, 0],
            'ind2': [0, 1],
            'ind3': [1, 1]
        }, index=['var1', 'var2'])
        
        expr = pd.DataFrame({
            'ind1': [5.5],
            'ind2': [6.2],
            'ind3': [4.8]
        })
        
        sad, c0 = nocovar_linear_estimate(h1, h2, expr)
        
        # Check that we get estimates for each variant
        assert 'var1' in sad
        assert 'var2' in sad
        assert isinstance(c0, (list, np.ndarray))


class TestIntegration:
    """Integration tests using real test data"""
    
    def test_full_pipeline_load(self):
        """Test loading all required data files"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=False)
        expr_file = read_expressions("tests/dataset/test_expressions.csv.gz", eqtl_file, args)
        
        expr_file.index = expr_file.Name
        expr_file = expr_file.drop(columns=["Name"]).astype(float)
        
        # Filter genes
        filtered_eqtl = drop_noexpr_genes(eqtl_file, expr_file)
        
        # Initialize output
        output_df = init_output_df(filtered_eqtl, args)
        
        assert not output_df.empty
        assert 'log2_aFC' in output_df.columns
        assert len(output_df) > 0
        
    def test_normalized_expressions_loading(self):
        """Test loading non-normalized expression data"""
        eqtl_file = read_eqtls("tests/dataset/test_eqtls.tsv")
        args = MockArgs(gct=True)
        
        expr_file = read_expressions("tests/dataset/test_expressions_nonorm.csv.gz", eqtl_file, args)
        
        assert not expr_file.empty
        assert "Name" in expr_file.columns


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
