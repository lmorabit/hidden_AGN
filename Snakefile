rule generateVmaxes:
	input:
                "src/static/cochrane_2023_tableA1.csv",
                "src/static/kondapally_2022_table1.csv",
		"src/static/lockman_03_matched_inMOC_inHR.fits",
                "src/static/lockman_DR1_rms_masked.fits",
		"src/static/en1_03_matched_inMOC_inHR.fits",
                "src/static/en1_DR1_rms_masked.fits"
	output:
		"src/data/lockman_vmaxes_zmin0.003_zmax0.3.fits",
		"src/data/en1_vmaxes_zmin0.003_zmax0.3.fits",
		"src/data/lockman_vmaxes_zmin0.5_zmax1.0.fits",
		"src/data/lockman_vmaxes_zmin1.0_zmax1.5.fits",
		"src/data/lockman_vmaxes_zmin1.5_zmax2.0.fits",
		"src/data/lockman_vmaxes_zmin2.0_zmax2.5.fits",
		"src/data/lockman_vmaxes_zmin0.1_zmax0.4.fits",
		"src/data/lockman_vmaxes_zmin0.4_zmax0.6.fits",
		"src/data/lockman_vmaxes_zmin0.6_zmax0.8.fits",
		"src/data/lockman_vmaxes_zmin0.8_zmax1.0.fits",
		"src/data/lockman_vmaxes_zmin1.0_zmax1.3.fits",
		"src/data/lockman_vmaxes_zmin1.3_zmax1.6.fits",
		"src/data/lockman_vmaxes_zmin1.6_zmax2.0.fits",
		"src/data/lockman_vmaxes_zmin2.5_zmax3.3.fits",
		"src/data/lockman_vmaxes_zmin3.3_zmax4.6.fits",
		"src/data/lockman_vmaxes_zmin4.6_zmax5.7.fits",
		"src/data/en1_vmaxes_zmin0.5_zmax1.0.fits",
		"src/data/en1_vmaxes_zmin1.0_zmax1.5.fits",
		"src/data/en1_vmaxes_zmin1.5_zmax2.0.fits",
		"src/data/en1_vmaxes_zmin2.0_zmax2.5.fits",
		"src/data/en1_vmaxes_zmin0.1_zmax0.4.fits",
		"src/data/en1_vmaxes_zmin0.4_zmax0.6.fits",
		"src/data/en1_vmaxes_zmin0.6_zmax0.8.fits",
		"src/data/en1_vmaxes_zmin0.8_zmax1.0.fits",
		"src/data/en1_vmaxes_zmin1.0_zmax1.3.fits",
		"src/data/en1_vmaxes_zmin1.3_zmax1.6.fits",
		"src/data/en1_vmaxes_zmin1.6_zmax2.0.fits",
		"src/data/en1_vmaxes_zmin2.5_zmax3.3.fits",
		"src/data/en1_vmaxes_zmin3.3_zmax4.6.fits",
		"src/data/en1_vmaxes_zmin4.6_zmax5.7.fits"
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
                "src/data/lockman_vmaxes_zmin0.003_zmax0.3.fits",
                "src/data/en1_vmaxes_zmin0.003_zmax0.3.fits"
        output:
                "deep_fields_RLFs.png"
        conda:
                "environment.yml"
        script:
                "src/scripts/comparison_plot.py"

