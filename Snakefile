rule generateVmaxes:
	input:
                "src/static/cochrane_2023_tableA1.csv",
                "src/static/kondapally_2022_table1.csv",
		"src/static/lockman_03_matched_inMOC_inHR.fits",
                "src/static/lockman_DR1_rms_masked.fits",
		"src/static/en1_03_matched_inMOC_inHR.fits",
                "src/static/en1_DR1_rms_masked.fits"
	output:
		"src/data/lockman_vmaxes.fits",
		"src/data/en1_vmaxes.fits"
	conda:
		"environment.yml"
	script:
		"src/scripts/generateVmaxes.py"
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
		"src/static/kondapally_2022_table2.csv",
		"src/static/cochrane_2023_table1.csv",
		"src/data/lockman_vmaxes.fits",
	        "src/data/en1_vmaxes.fits"
        output:
                "mauch_sadler_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/comparison_plot.py"

