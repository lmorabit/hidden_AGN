rule generate_vmax:
	input:
                "src/data/kondapally_2022_table1.csv",
                "src/data/cochrane_2023_tableA1.csv",
		"src/static/en1_03_matched_inMOC_inHR.fits",
		"src/static/lockman_03_matched_inMOC_inHR.fits",
		"src/static/en1_DR1_rms_masked.fits",
		"src/static/lockman_DR1_rms_masked.fits"
	output:
		directory("src/data/vmaxes")
	cache:
		True
	conda:
		"environment.yml"
	script:
		"src/scripts/generateVmaxes.py"
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
rule rlf_evolution:
	input:
                "src/static/redshift_bins.csv",
		"src/data/lockman_vmaxes_zmin0.003_zmax0.3.fits",
                "src/data/en1_vmaxes_zmin0.003_zmax0.3.fits",
                "src/data/lockman_vmaxes_zmin0.5_zmax1.0.fits",
                "src/data/lockman_vmaxes_zmin1.0_zmax1.5.fits",
                "src/data/lockman_vmaxes_zmin1.5_zmax2.0.fits",
                "src/data/lockman_vmaxes_zmin2.0_zmax2.5.fits",
                "src/data/en1_vmaxes_zmin0.5_zmax1.0.fits",
                "src/data/en1_vmaxes_zmin1.0_zmax1.5.fits",
                "src/data/en1_vmaxes_zmin1.5_zmax2.0.fits",
                "src/data/en1_vmaxes_zmin2.0_zmax2.5.fits"
	output:
		"RLF_evolution.png",
                "src/output/integrated_differences.txt"
	conda:
		"environment.yml"
	script:
		"src/scripts/RLF_evolution.py"
rule calculate_vars:
	input: 
		"src/static/en1_03_matched_inMOC_inHR.fits",
		"src/static/lockman_03_matched_inMOC_inHR.fits",
	output:
                "src/tex/output/en1_detectable.txt",
                "src/tex/output/flowchart_numbers.txt",
                "src/tex/output/lockman_detectable.txt"
	conda:
		"environment.yml"
	script:
		"src/scripts/flowchart_numbers.py"
