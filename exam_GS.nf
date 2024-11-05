nextflow.enable.dsl = 2

params.accession = "M21012"
params.out = "${launchDir}/output"
params.storeDir="${launchDir}/cache"


process downloadAccession {
	storeDir params.storeDir
	input:
		val accession
	output:
		path "${params.accession}.fasta"
		
	script:
	"""
	wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${params.accession}&rettype=fasta&retmode=text" -O ${params.accession}.fasta
	"""
}


process downloadSequences {
	storeDir params.storeDir
	output:
		path "sequences.fasta"
		
	script:
	"""
	wget https://gitlab.com/dabrowskiw/cq-examples/-/raw/master/data/hepatitis_combined.fasta?inline=false -O sequences.fasta
	"""
}


process combineFasta {
	publishDir params.out, mode: "copy", overwrite: true
	input:
		path infile
	output:
		path "${params.accession}_combined.fasta"
	
	script:
	"""
	cat *.fasta > ${params.accession}_combined.fasta
	"""
}


process mafft {
	publishDir params.out, mode: "copy", overwrite: true
    container "https://depot.galaxyproject.org/singularity/mafft%3A7.525--h031d066_1"
    input:
        path infile
    output:
        path "${params.accession}_mafft.fasta"
		
    script:
    """
    mafft --auto ${infile} > ${params.accession}_mafft.fasta
    """
}


process trimAL {
	publishDir params.out, mode: "copy", overwrite: true
	container "https://depot.galaxyproject.org/singularity/trimal%3A1.5.0--h4ac6f70_1"
	input:
		path infile
	output:
		path "${infile}*"
		
	script:
	"""
	trimal -in $infile -out ${infile}.trimmed.fasta -htmlout ${infile}_report.html -automated1 
	"""
}


workflow {
	accession_channel = Channel.from(params.accession)
    download_accession_channel = downloadAccession(accession_channel)
	
	download_sequences_channel = downloadSequences()
	
	concat_channel = download_accession_channel.concat(download_sequences_channel).collect()
	combine_fasta_channel = combineFasta(concat_channel)
	alignment_channel = mafft(combine_fasta_channel)
	trim_channel = trimAL(alignment_channel)
}