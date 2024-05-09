rule generateVmaxes:
	input:
		"src/static/lockman_03_matched_inMOC_inHR.fits",
                "src/static/lockman_DR1_rms_masked.fits",
                "src/static/cochrane_2023_tableA1.csv",
                "src/static/kondapally_2022_table1.csv"
	output:
		"src/data/lockman_6arcsec_vmaxes.fits"
	conda:
		"environment.yml"
	script:
		"src/scripts/generateVmaxes.py"

rule testplot:
        input:
                "src/static/RLFS_50MYR_SF.csv",
                "src/static/RLFS_50MYR_AGN.csv"
        output:
                "test.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/testplot.py"
rule completeness:
	input:
		"src/static/cochrane_2023_table1.csv",
		"src/static/cochrane_2023_tableA1.csv",
		"src/static/kondapally_2022_table1.csv",
		"src/static/kondapally_2022_table2.csv"
	output:
		"completeness.png"
	conda:
		"environment.yml"
	script:
		"src/scripts/completeness.py"
rule comparison_plot:
        input:
                "src/data/lockman_6arcsec_vmaxes.fits",
		"src/static/kondapally_2022_table2.csv",
		"src/static/cochrane_2023_table1.csv"
        output:
                "mauch_sadler_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/comparison_plot.py"

