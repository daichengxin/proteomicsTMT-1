#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/proteomicstmt
========================================================================================
 nf-core/proteomicstmt Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/proteomicstmt
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/proteomicstmt --spectra '*.mzML' --database '*.fasta' -profile docker

    Main arguments:
      --input                       Path/URI to PRIDE Sample to data relation format file (SDRF) OR path to input spectra as mzML or Thermo Raw

      For SDRF:
      --root_folder                 (Optional) If given, looks for the filenames in the SDRF in this folder, locally
      --local_input_type            (Optional) If given and 'root_folder' was specified, it overwrites the filetype in the SDRF for local lookup and matches only the basename.

      For mzML/raw files:
      --expdesign                   (Optional) Path to an experimental design file (if not given, it assumes unfractionated, unrelated samples)

      And:
      --database                    Path to input protein database as fasta

    Decoy database:
      --add_decoys                  Add decoys to the given fasta
      --decoy_affix                 The decoy prefix or suffix used or to be used (default: DECOY_)
      --affix_type                  Prefix (default) or suffix (WARNING: Percolator only supports prefices)

    Isobaric analyze:
      --label   					label method (TMT6plex,TMT10plex,TMT11plex)
      --fragment_method				dissolve method (HCD,CID)
      --min_precursor_intensity     Minimum intensity of the precursor to be extracted. (Default 1.0)
      --normalization 				Enable normalization of channel intensities with respect to the reference channel.(default false)
      --reference_channel			Number of the reference channel
      --isotope_correction			Enable isotope correction (highly recommended,default true)


    Database Search:
      --search_engines               Which search engine: "comet" (default) or "msgf"
      --enzyme                      Enzymatic cleavage (e.g. 'unspecific cleavage' or 'Trypsin' [default], see OpenMS enzymes)
      --num_enzyme_termini          Specify the termini where the cleavage rule has to match (default:
                                         'fully' valid: 'semi', 'fully')
      --num_hits                    Number of peptide hits per spectrum (PSMs) in output file (default: '1')
      --fixed_mods                  Fixed modifications ('Carbamidomethyl (C)', see OpenMS modifications)
      --variable_mods               Variable modifications ('Oxidation (M)', see OpenMS modifications)
      --precursor_mass_tolerance    Mass tolerance of precursor mass (default: 5)
      --precursor_mass_tolerance_unit Da or ppm (default: ppm)
      --fragment_mass_tolerance     Mass tolerance for fragment masses (currently only controls Comets fragment_bin_tol) (default: 0.03)
      --fragment_mass_tolerance_unit Da or ppm (default: Da)
      --allowed_missed_cleavages    Allowed missed cleavages (default: 2)
      --min_precursor_charge        Minimum precursor ion charge (default: 2)
      --max_precursor_charge        Maximum precursor ion charge (default: 4)
      --min_peptide_length          Minimum peptide length to consider (default: 6)
      --max_peptide_length          Maximum peptide length to consider (default: 40)
      --instrument                  Type of instrument that generated the data (currently only 'high_res' [default] and 'low_res' supported)
      --protocol                    Used labeling or enrichment protocol (if any, phospho,iTRAQ,iTRAQ_phospho,TMT)
      --fragment_method             Used fragmentation method (currently unused since we let the search engines consider all MS2 spectra and let them determine from the spectrum metadata)
      --max_mods                    Maximum number of modifications per peptide. If this value is large, the search may take very long
      --db_debug                    Debug level during database search

      //TODO probably also still some options missing. Try to consolidate them whenever the two search engines share them

    Peak picking:
      --openms_peakpicking          Use the OpenMS PeakPicker to ADDITIONALLY pick the spectra before the search. This is usually done
                                    during conversion already. Only activate if something goes wrong.
      --peakpicking_inmemory        Perform OpenMS peakpicking in-memory. Needs at least the size of the mzML file as RAM but is faster. default: false
      --peakpicking_ms_levels       Which MS levels to pick. default: [] which means auto-convert all non-centroided

    Peptide Re-indexing:
      --IL_equivalent               Should isoleucine and leucine be treated interchangeably? Default: true
      --allow_unmatched             Ignore unmatched peptides (Default: false; only activate if you double-checked all other settings)

    PSM Rescoring:
      --posterior_probabilities     How to calculate posterior probabilities for PSMs:
                                    "percolator" = Re-score based on PSM-feature-based SVM and transform distance
                                        to hyperplane for posteriors
                                    "fit_distributions" = Fit positive and negative distributions to scores
                                        (similar to PeptideProphet)
      --rescoring_debug             Debug level during PSM rescoring
      --psm_pep_fdr_cutoff          FDR cutoff on PSM level (or potential peptide level; see Percolator options) before going into
                                    feature finding, map alignment and inference.

      Percolator specific:
      --train_FDR                   False discovery rate threshold to define positive examples in training. Set to testFDR if 0
      --test_FDR                    False discovery rate threshold for evaluating best cross validation result and reported end result
      --percolator_fdr_level        Level of FDR calculation ('peptide-level-fdrs' or 'psm-level-fdrs')
      --description_correct_features Description of correct features for Percolator (0, 1, 2, 4, 8, see Percolator retention time and calibration)
      --generic_feature_set         Use only generic (i.e. not search engine specific) features. Generating search engine specific
                                    features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.
      --subset_max_train            Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other
                                    PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.
      --klammer                     Retention time features are calculated as in Klammer et al. instead of with Elude

      Distribution specific:
      --outlier_handling            How to handle outliers during fitting:
                                    - ignore_iqr_outliers (default): ignore outliers outside of 3*IQR from Q1/Q3 for fitting
                                    - set_iqr_to_closest_valid: set IQR-based outliers to the last valid value for fitting
                                    - ignore_extreme_percentiles: ignore everything outside 99th and 1st percentile (also removes equal values like potential censored max values in XTandem)
                                    - none: do nothing
      --top_hits_only               Use only the top hits for fitting

      //TODO add more options for rescoring part


    ConsensusID:
      --consensusid_algorithm       Choose method to combine probabilities from multiple search engines (if used). Valid: best, worst, average, rank, PEPMatrix, PEPIons (Default: best)
      --min_consensus_support       Choose ratio of ADDITIONAL evidence for a peptide ID of a spectrum. Varies across methods. See documentation for further info. (Default: 0)
      --consensusid_considered_top_hits Number of top hits per spectrum considered for consensus scoring. (Default: 0 = all)


    IDFilter:
      --protein_level_fdr_cutoff      Protein level FDR cutoff (this affects and chooses the peptides used for quantification)


    IDMapper:
      --rt_tolerance				RT tolerance (in seconds) for the matching of peptide identifications and consensus features (Default 5.0).
      --mz_tolerance				m/z tolerance (in ppm or Da) for the matching of peptide identifications and consensus features(Default 20.0).
      --mz_measure					Unit of 'mz_tolerance'. ('ppm', 'Da',Default 'ppm')
      --mz_reference				Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used.
									If 'peptide',masses are computed from the sequences of peptide hits.(Defalut 'peptide')

	  FileMerger:
	    --annotate_file_origin		Store the original filename in each feature using meta value "file_origin".(Default false).



    Inference:
      --greedy_group_resolution		Default none.
      --top_PSMs					Consider only top X PSMs per spectrum. 0 considers all.(Default 1).

    IDConflictResolver:
      --resolve_between_features	A map may contain multiple features with both identical (possibly modified i.e. not stripped) sequence and charge state. 							   The feature with the 'highest intensity' is very likely the most reliable one.Default(off),highest_intensity

    ProteinQuantifier:
      --top 						Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all, Default 3)
      --average						Averaging method used to compute protein abundances from peptide abundances.
      								(median,mean,weighted_mean,sum Default median).
      --best_charge_and_fraction	Distinguish between fraction and charge states of a peptide.(Default false).
      --ratios						Add the log2 ratios of the abundance values to the output.(Default false)
      --normalize					Scale peptide abundances so that medians of all samples are equal.(Default false)
      --fix_peptides				Use the same peptides for protein quantification across all samples.(Default false)



    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name

if (!(workflow.runName == ~/[a-z]+_[a-z]+/)) {
	custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}


// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// Validate input
if (isCollectionOrArray(params.input))
{
  tocheck = params.input[0]
} else {
  tocheck = params.input
}


sdrf_file = null

if (tocheck.toLowerCase().endsWith("sdrf") || tocheck.toLowerCase().endsWith("tsv")) {
  sdrf_file = params.input
} else if (tocheck.toLowerCase().endsWith("mzml") || tocheck.toLowerCase().endsWith("raw")) {
  spectra_files = params.input
} else {
	log.error "EITHER spectra data (mzml/raw) OR an SDRF needs to be  provided as input."; exit 1
}

params.database = params.database ?: { log.error "No protein database provided. Make sure you have used the '--database' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Create a channel for input files
 */

 //Filename        FixedModifications      VariableModifications   Label   PrecursorMassTolerance  PrecursorMassToleranceUnit      FragmentMassTolerance   DissociationMethod      Enzyme


if (!sdrf_file) {
  ch_spectra = Channel.fromPath(spectra_files, checkIfExists: true)
  ch_spectra
  .multiMap{ it -> id = it.toString().md5()
  					isobaricanalyzer_settings: tuple(id,
  									params.label,
  									params.fragment_method)
                    comet_settings: msgf_settings: tuple(id,
                                    params.fixed_mods,
                                    params.variable_mods,
                                    "", //labelling modifications currently not supported
                                    params.precursor_mass_tolerance,
                                    params.precursor_mass_tolerance_unit,
                                    params.fragment_mass_tolerance,
                                    params.fragment_mass_tolerance_unit,
                                    params.fragment_method,
                                    params.enzyme)
                    idx_settings: tuple(id,
                                    params.enzyme)
                    luciphor_settings:
                                  tuple(id,
                                    params.fragment_method)
                    mzmls: tuple(id,it)}
  .set{ch_sdrf_config}
} else {
  ch_sdrf = Channel.fromPath(sdrf_file, checkIfExists: true)
  /*
   * STEP 0 - SDRF parsing
   */
  process sdrf_parsing {

      publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

      input:
       file sdrf from ch_sdrf

      output:
       file "experimental_design.tsv" into ch_expdesign, ch_pro_quant_exp
       file "openms.tsv" into ch_sdrf_config_file

      when:
        sdrf_file

      script:
       """
       ## -t2 since the one-table format parser is broken in OpenMS2.5
       ## -l for legacy behavior to always add sample columns
       parse_sdrf convert-openms -t2 -l -s ${sdrf} > sdrf_parsing.log
       """
  }

  //TODO use header and reference by col name instead of index
  ch_sdrf_config_file
  .splitCsv(skip: 1, sep: '\t')
  .multiMap{ row -> id = row.toString().md5()
  					isobaricanalyzer_settings: tuple(id,
  									row[4],
  									row[9])
                    comet_settings: msgf_settings: tuple(id,
                                    row[2],
                                    row[3],
                                    row[4],
                                    row[5],
                                    row[6],
                                    row[7],
                                    row[8],
                                    row[9],
                                    row[10])
                    idx_settings: tuple(id,
                                    row[10])
                    luciphor_settings:
                                  tuple(id,
                                    row[9])
                    mzmls: tuple(id, !params.root_folder ?
                                    row[0] :
                                    params.root_folder + "/" + (params.local_input_type ?
                                        row[1].take(row[1].lastIndexOf('.')) + '.' + params.local_input_type :
                                        row[1]))}
  .set{ch_sdrf_config}
}



ch_db_for_decoy_creation = Channel.fromPath(params.database)

// overwrite experimental design if given additionally to SDRF
//TODO think about that
if (params.expdesign)
{
    Channel
        .fromPath(params.expdesign)
        .set { ch_pro_quant_exp }
}


ch_sdrf_config.mzmls
.branch {
        raw: hasExtension(it[1], 'raw')
        mzML: hasExtension(it[1], 'mzML')
}
.set {branched_input}



//TODO we could also check for outdated mzML versions and try to update them
branched_input.mzML
.branch {
    nonIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return false;
                        }
                        return true;
                    }
    inputIndexedMzML: file(it[1]).withReader {
                        f = it;
                        1.upto(5) {
                            if (f.readLine().contains("indexedmzML")) return true;
                        }
                        return false;
                    }
}
.set {branched_input_mzMLs}

//Push raw files through process that does the conversion, everything else directly to downstream Channel with mzMLs

//This piece only runs on data that is a.) raw and b.) needs conversion
//mzML files will be mixed after this step to provide output for downstream processing - allowing you to even specify mzMLs and RAW files in a mixed mode as input :-)

/*
 * STEP 0.1 - Raw file conversion
 */
process raw_file_conversion {

    label 'process_low'
    label 'process_single_thread'


    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(rawfile) from branched_input.raw

    output:
     tuple mzml_id, file("*.mzML") into mzmls_converted

    script:
     """
     ThermoRawFileParser.sh -i=${rawfile} -f=2 -o=./ > ${rawfile}_conversion.log
     """
}

/*
 * STEP 0.2 - MzML indexing
 */
process mzml_indexing {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(mzmlfile) from branched_input_mzMLs.nonIndexedMzML

    output:
     tuple mzml_id, file("out/*.mzML") into mzmls_indexed
     file "*.log"

    script:
     """
     mkdir out
     FileConverter -in ${mzmlfile} -out out/${mzmlfile.baseName}.mzML > ${mzmlfile.baseName}_mzmlindexing.log
     """
}

//Mix the converted raw data with the already supplied mzMLs and push these to the same channels as before

if (params.openms_peakpicking)
{
  branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).set{mzmls_pp}
  (mzmls_isobaric_analyzer, mzmls_msgf, mzmls_comet, mzmls_luciphor) = [Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty()]
}
else
{
  branched_input_mzMLs.inputIndexedMzML.mix(mzmls_converted).mix(mzmls_indexed).into{mzmls_isobaric_analyzer; mzmls_comet; mzmls_msgf; mzmls_luciphor}
  mzmls_pp = Channel.empty()
}

//Fill the channels with empty Channels in case that we want to add decoys. Otherwise fill with output from database.
(searchengine_in_db_msgf, searchengine_in_db_comet, pepidx_in_db, plfq_in_db) = ( params.add_decoys
                    ? [ Channel.empty(), Channel.empty(), Channel.empty(), Channel.empty() ]
                    : [ Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database), Channel.fromPath(params.database) ] )

//Add decoys if params.add_decoys is set appropriately
process generate_decoy_database {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     file(mydatabase) from ch_db_for_decoy_creation

    output:
     file "${mydatabase.baseName}_decoy.fasta" into searchengine_in_db_decoy_msgf, searchengine_in_db_decoy_comet, pepidx_in_db_decoy, plfq_in_db_decoy
     file "*.log"

    when:
     params.add_decoys

    script:
     """
     DecoyDatabase  -in ${mydatabase} \\
                 -out ${mydatabase.baseName}_decoy.fasta \\
                 -decoy_string ${params.decoy_affix} \\
                 -decoy_string_position ${params.affix_type} \\
                 > ${mydatabase.baseName}_decoy_database.log
     """
}


process openms_peakpicker {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, path(mzml_file) from mzmls_pp

    when:
      params.openms_peakpicking

    output:
     tuple mzml_id, file("out/${mzml_file.baseName}.mzML") into isobaric_analyzer_picked, mzmls_msgf_picked, mzmls_comet_picked, mzmls_plfq_picked
     file "*.log"

    script:
     in_mem = params.peakpicking_inmemory ? "inmemory" : "lowmemory"
     lvls = params.peakpicking_ms_levels ? "-algorithm:ms_levels ${params.peakpicking_ms_levels}" : ""
     """
     mkdir out
     PeakPickerHiRes -in ${mzml_file} \\
                     -out out/${mzml_file.baseName}.mzML \\
                     -threads ${task.cpus} \\
                     -debug ${params.pp_debug} \\
                     -processOption ${in_mem} \\
                     ${lvls} \\
                     > ${mzml_file.baseName}_pp.log
     """
}

if (params.enzyme == "unspecific cleavage")
{
  params.num_enzyme_termini == "none"
}

pepidx_num_enzyme_termini = params.num_enzyme_termini
if (params.num_enzyme_termini == "fully")
{
  pepidx_num_enzyme_termini = "full"
}


process isobaric_analyzer {
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	  tuple mzml_id, path(mzml_file), label,  diss_meth from mzmls_isobaric_analyzer.mix(isobaric_analyzer_picked).join(ch_sdrf_config.isobaricanalyzer_settings)

	when:
	  params.isobaric_analyzer


	output:
	  tuple mzml_id, file("${mzml_file.baseName}_iso.consensusXML") into id_files_consensusXML
      file "*.log"

	script:
	  if (diss_meth == 'HCD') diss_meth = 'High-energy collision-induced dissociation'
	  else if (diss_meth == 'CID') diss_meth = 'Collision-induced dissociation'
	  else if (diss_meth == 'ETD') diss_meth = 'Electron transfer dissociation'
	  else if (diss_meth == 'ECD') diss_meth = 'Electron capture dissociation'
	  """
	  IsobaricAnalyzer -type ${label} \\
	  				   -in ${mzml_file} \\
	  				   -threads ${task.cpus} \\
	  				   -quantification:normalization \\
	  				   -extraction:select_activation "${diss_meth}" \\
	  				   -extraction:min_reporter_intensity ${params.min_reporter_intensity} \\
	  				   -extraction:min_precursor_purity ${params.min_precursor_purity} \\
	  				   -extraction:precursor_isotope_deviation ${params.precursor_isotope_deviation} \\
	  				   -debug ${params.iso_debug} \\
	  				   -out ${mzml_file.baseName}_iso.consensusXML \\
	  				   > ${mzml_file.baseName}_isob.log
	  """


}


process search_engine_msgf {

    label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------- WARNING: THIS IS A HACK. IT JUST DOES NOT WORK IF THIS PROCESS IS RETRIED -----------------------------
    // ---------------------------------------------------------------------------------------------------------------------
    // I actually dont know, where else this would be needed.
    errorStrategy 'terminate'

    input:
     tuple file(database), mzml_id, path(mzml_file), fixed, variable, label, prec_tol, prec_tol_unit, frag_tol, frag_tol_unit, diss_meth, enzyme from searchengine_in_db_msgf.mix(searchengine_in_db_decoy_msgf).combine(mzmls_msgf.mix(mzmls_msgf_picked).join(ch_sdrf_config.msgf_settings))

     // This was another way of handling the combination
     //file database from searchengine_in_db.mix(searchengine_in_db_decoy)
     //each file(mzml_file) from mzmls
    when:
      params.search_engines.contains("msgf")

    output:
     tuple mzml_id, file("${mzml_file.baseName}_msgf.idXML") into id_files_msgf
     file "*.log"

    script:
      if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
      else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
      else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
      else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
      else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'

      if ((frag_tol.toDouble() < 50 && frag_tol_unit == "ppm") || (frag_tol.toDouble() < 0.1 && frag_tol_unit == "Da"))
      {
        inst = params.instrument ?: "high_res"
      } else {
        inst = params.instrument ?: "low_res"
      }
     """
     MSGFPlusAdapter -in ${mzml_file} \\
                     -out ${mzml_file.baseName}_msgf.idXML \\
                     -threads ${task.cpus} \\
                     -java_memory ${task.memory.toMega()} \\
                     -database "${database}" \\
                     -instrument ${inst} \\
                     -protocol "${params.protocol}" \\
                     -matches_per_spec ${params.num_hits} \\
                     -min_precursor_charge ${params.min_precursor_charge} \\
                     -max_precursor_charge ${params.max_precursor_charge} \\
                     -min_peptide_length ${params.min_peptide_length} \\
                     -max_peptide_length ${params.max_peptide_length} \\
                     -enzyme "${enzyme}" \\
                     -tryptic ${params.num_enzyme_termini} \\
                     -precursor_mass_tolerance ${prec_tol} \\
                     -precursor_error_units ${prec_tol_unit} \\
                     -fixed_modifications ${fixed.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -variable_modifications ${variable.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                     -max_mods ${params.max_mods} \\
                     -debug ${params.db_debug} \\
                     > ${mzml_file.baseName}_msgf.log
     """
}


process search_engine_comet {

    label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    // ---------------------------------------------------------------------------------------------------------------------
    // ------------- WARNING: THIS IS A HACK. IT JUST DOES NOT WORK IF THIS PROCESS IS RETRIED -----------------------------
    // ---------------------------------------------------------------------------------------------------------------------
    // I actually dont know, where else this would be needed.
    errorStrategy 'terminate'
    input:
     tuple file(database), mzml_id, path(mzml_file), fixed, variable, label, prec_tol, prec_tol_unit, frag_tol, frag_tol_unit, diss_meth, enzyme from searchengine_in_db_comet.mix(searchengine_in_db_decoy_comet).combine(mzmls_comet.mix(mzmls_comet_picked).join(ch_sdrf_config.comet_settings))

    when:
      params.search_engines.contains("comet")

    output:
     tuple mzml_id, file("${mzml_file.baseName}_comet.idXML") into id_files_comet
     file "*.log"

    //TODO we currently ignore the activation_method param to leave the default "ALL" for max. compatibility
    script:
     if (frag_tol_unit == "ppm") {
       // Note: This uses an arbitrary rule to decide if it was hi-res or low-res
       // and uses Comet's defaults for bin size, in case unsupported unit "ppm" was given.
       if (frag_tol.toDouble() < 50) {
         bin_tol = "0.015"
         bin_offset = "0.0"
         inst = params.instrument ?: "high_res"
       } else {
         bin_tol = "0.50025"
         bin_offset = "0.4"
         inst = params.instrument ?: "low_res"
       }
       log.warn "The chosen search engine Comet does not support ppm fragment tolerances. We guessed a " + inst +
         " instrument and set the fragment_bin_tolerance to " + bin_tol
     } else {
       bin_tol = frag_tol.toDouble() / 2.0
       bin_offset = frag_tol.toDouble() < 0.1 ? "0.0" : "0.4"
       if (!params.instrument)
       {
         inst = frag_tol.toDouble() < 0.1 ? "high_res" : "low_res"
       } else {
         inst = params.instrument
       }
     }

     // for consensusID the cutting rules need to be the same. So we adapt to the loosest rules from MSGF
     // TODO find another solution. In ProteomicsLFQ we re-run PeptideIndexer (remove??) and if we
     // e.g. add XTandem, after running ConsensusID it will lose the auto-detection ability for the
     // XTandem specific rules.
     if (params.search_engines.contains("msgf"))
     {
        if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
        else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
        else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
        else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
        else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'
     }
     """
     CometAdapter  -in ${mzml_file} \\
                   -out ${mzml_file.baseName}_comet.idXML \\
                   -threads ${task.cpus} \\
                   -database "${database}" \\
                   -instrument ${inst} \\
                   -missed_cleavages ${params.allowed_missed_cleavages} \\
                   -num_hits ${params.num_hits} \\
                   -num_enzyme_termini ${params.num_enzyme_termini} \\
                   -enzyme "${enzyme}" \\
                   -precursor_charge ${params.min_precursor_charge}:${params.max_precursor_charge} \\
                   -fixed_modifications ${fixed.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                   -variable_modifications ${variable.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                   -max_variable_mods_in_peptide ${params.max_mods} \\
                   -precursor_mass_tolerance ${prec_tol} \\
                   -precursor_error_units ${prec_tol_unit} \\
                   -fragment_mass_tolerance ${bin_tol} \\
                   -fragment_bin_offset ${bin_offset} \\
                   -debug ${params.db_debug} \\
                   > ${mzml_file.baseName}_comet.log
     """
}


process index_peptides {

    label 'process_low'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file), val(enzyme), file(database) from id_files_msgf.mix(id_files_comet).combine(ch_sdrf_config.idx_settings, by: 0).combine(pepidx_in_db.mix(pepidx_in_db_decoy))

    output:
     tuple mzml_id, file("${id_file.baseName}_idx.idXML") into id_files_idx_ForPerc, id_files_idx_ForIDPEP, id_files_idx_ForIDPEP_noFDR
     file "*.log"

    script:
     def il = params.IL_equivalent ? '-IL_equivalent' : ''
     def allow_um = params.allow_unmatched ? '-allow_unmatched' : ''
     // see comment in CometAdapter. Alternative here in PeptideIndexer is to let it auto-detect the enzyme by not specifying.
     if (params.search_engines.contains("msgf"))
     {
        if (enzyme == 'Trypsin') enzyme = 'Trypsin/P'
        else if (enzyme == 'Arg-C') enzyme = 'Arg-C/P'
        else if (enzyme == 'Asp-N') enzyme = 'Asp-N/B'
        else if (enzyme == 'Chymotrypsin') enzyme = 'Chymotrypsin/P'
        else if (enzyme == 'Lys-C') enzyme = 'Lys-C/P'
     }
     """
     PeptideIndexer -in ${id_file} \\
                    -out ${id_file.baseName}_idx.idXML \\
                    -threads ${task.cpus} \\
                    -fasta ${database} \\
                    -enzyme:name "${enzyme}" \\
                    -enzyme:specificity ${pepidx_num_enzyme_termini} \\
                    ${il} \\
                    ${allow_um} \\
                    > ${id_file.baseName}_index_peptides.log
     """
}


// ---------------------------------------------------------------------
// Branch a) Q-values and PEP from Percolator

process extract_percolator_features {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_ForPerc

    output:
     tuple mzml_id, file("${id_file.baseName}_feat.idXML") into id_files_idx_feat
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    script:
     """
     PSMFeatureExtractor -in ${id_file} \\
                         -out ${id_file.baseName}_feat.idXML \\
                         -threads ${task.cpus} \\
                         > ${id_file.baseName}_extract_percolator_features.log
     """
}


//Note: from here, we do not need any settings anymore. so we can skip adding the mzml_id to the channels
//TODO find a way to run across all runs merged
process percolator {

    //TODO Actually it heavily depends on the subset_max_train option and the number of IDs
    // would be cool to get an estimate by parsing the number of IDs from previous tools.
    label 'process_medium'
    //Since percolator 3.5 it allows for 27 parallel tasks
    cpus { check_max( 27, 'cpus' ) }

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/raw_ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_feat

    output:
     tuple mzml_id, file("${id_file.baseName}_perc.idXML"), val("MS:1001491") into id_files_perc, id_files_perc_consID
     file "*.log"

    when:
     params.posterior_probabilities == "percolator"

    // NICE-TO-HAVE: the decoy-pattern is automatically detected from PeptideIndexer.
    // Parse its output and put the correct one here.
    script:
      if (params.klammer && params.description_correct_features == 0) {
          log.warn('Klammer was specified, but description of correct features was still 0. Please provide a description of correct features greater than 0.')
          log.warn('Klammer will be implicitly off!')
      }

      // currently post-processing-tdc is always set since we do not support separate TD databases
      """
      ## Percolator does not have a threads parameter. Set it via OpenMP env variable,
      ## to honor threads on clusters
      OMP_NUM_THREADS=${task.cpus} PercolatorAdapter \\
                          -in ${id_file} \\
                          -out ${id_file.baseName}_perc.idXML \\
                          -threads ${task.cpus} \\
                          -subset_max_train ${params.subset_max_train} \\
                          -decoy_pattern ${params.decoy_affix} \\
                          -post_processing_tdc \\
                          -score_type pep \\
                          > ${id_file.baseName}_percolator.log
      """
}


// ---------------------------------------------------------------------
// Branch b) Q-values and PEP from OpenMS

if(params.posterior_probabilities != "percolator" && params.search_engines.split(",").size() == 1)
{
  id_files_idx_ForIDPEP_noFDR = Channel.empty()
}
process fdr_idpep {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_ForIDPEP

    output:
     tuple mzml_id, file("${id_file.baseName}_fdr.idXML") into id_files_idx_ForIDPEP_FDR
     file "*.log"

    when:
     params.posterior_probabilities != "percolator" && params.search_engines.split(",").size() == 1

    script:
     """
     FalseDiscoveryRate -in ${id_file} \\
                        -out ${id_file.baseName}_fdr.idXML \\
                        -threads ${task.cpus} \\
                        -protein false \\
                        -algorithm:add_decoy_peptides \\
                        -algorithm:add_decoy_proteins \\
                        > ${id_file.baseName}_fdr.log
     """
}


//idpep picks the best scores for each search engine automatically. No switching needed after FDR.
process idpep {

    label 'process_low'
    // I think Eigen optimization is multi-threaded, so leave threads open

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/raw_ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file) from id_files_idx_ForIDPEP_FDR.mix(id_files_idx_ForIDPEP_noFDR)

    output:
     tuple mzml_id, file("${id_file.baseName}_idpep.idXML"), val("q-value_score") into id_files_idpep, id_files_idpep_consID
     file "*.log"

    when:
     params.posterior_probabilities != "percolator"

    script:
     """
     IDPosteriorErrorProbability    -in ${id_file} \\
                                    -out ${id_file.baseName}_idpep.idXML \\
                                    -fit_algorithm:outlier_handling ${params.outlier_handling} \\
                                    -threads ${task.cpus} \\
                                    > ${id_file.baseName}_idpep.log
     """
}


// ---------------------------------------------------------------------
// Main Branch

process idscoreswitcher_to_qval {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file), val(qval_score) from id_files_idpep.mix(id_files_perc)

    output:
     tuple mzml_id, file("${id_file.baseName}_switched.idXML") into id_files_noConsID_qval
     file "*.log"

    when:
     params.search_engines.split(",").size() == 1

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_switched.idXML \\
                        -threads ${task.cpus} \\
                        -old_score "Posterior Error Probability" \\
                        -new_score ${qval_score} \\
                        -new_score_type q-value \\
                        -new_score_orientation lower_better \\
                        > ${id_file.baseName}_scoreswitcher_qval.log
     """
}


process consensusid {

    label 'process_medium'
    //TODO could be easily parallelized
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/consensus_ids", mode: 'copy', pattern: '*.idXML'

    // we can drop qval_score in this branch since we have to recalculate FDR anyway
    input:
     tuple mzml_id, file(id_files_from_ses), val(qval_score) from id_files_idpep_consID.mix(id_files_perc_consID).groupTuple(size: params.search_engines.split(",").size())

    output:
     tuple mzml_id, file("${mzml_id}_consensus.idXML") into consensusids
     file "*.log"

    when:
     params.search_engines.split(",").size() > 1

    script:
     """
     ConsensusID -in ${id_files_from_ses} \\
                        -out ${mzml_id}_consensus.idXML \\
                        -per_spectrum \\
                        -threads ${task.cpus} \\
                        -algorithm ${params.consensusid_algorithm} \\
                        -filter:min_support ${params.min_consensus_support} \\
                        -filter:considered_hits ${params.consensusid_considered_top_hits} \\
                        > ${mzml_id}_consensusID.log
     """

}


process fdr_consensusid {

    label 'process_medium'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file) from consensusids

    output:
     tuple mzml_id, file("${id_file.baseName}_fdr.idXML") into consensusids_fdr
     file "*.log"

    when:
     params.search_engines.split(",").size() > 1

    script:
     """
     FalseDiscoveryRate -in ${id_file} \\
                        -out ${id_file.baseName}_fdr.idXML \\
                        -threads ${task.cpus} \\
                        -protein false \\
                        -algorithm:add_decoy_peptides \\
                        -algorithm:add_decoy_proteins \\
                        > ${id_file.baseName}_fdr.log
     """

}


process idfilter {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ids", mode: 'copy', pattern: '*.idXML'

    input:
     tuple mzml_id, file(id_file) from id_files_noConsID_qval.mix(consensusids_fdr)

    output:
     tuple mzml_id, file("${id_file.baseName}_filter.idXML") into id_filtered, id_filtered_luciphor
     file "*.log"

    script:
     """
     IDFilter -in ${id_file} \\
              -out ${id_file.baseName}_filter.idXML \\
              -threads ${task.cpus} \\
              -score:pep ${params.psm_pep_fdr_cutoff} \\
              > ${id_file.baseName}_idfilter.log
     """
}


ptmt_in_id = params.enable_mod_localization
                    ? Channel.empty()
                    : id_filtered

// TODO make luciphor pick its own score so we can skip this step
process idscoreswitcher_for_luciphor {

    label 'process_very_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(id_file) from id_filtered_luciphor

    output:
     tuple mzml_id, file("${id_file.baseName}_pep.idXML") into id_filtered_luciphor_pep
     file "*.log"

    when:
     params.enable_mod_localization

    script:
     """
     IDScoreSwitcher    -in ${id_file} \\
                        -out ${id_file.baseName}_pep.idXML \\
                        -threads ${task.cpus} \\
                        -old_score "q-value" \\
                        -new_score "Posterior Error Probability_score" \\
                        -new_score_type "Posterior Error Probability" \\
                        -new_score_orientation lower_better \\
                        > ${id_file.baseName}_switch_pep_for_luciphor.log
     """
}


process luciphor {

    label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
     tuple mzml_id, file(mzml_file), file(id_file), frag_method from mzmls_luciphor.join(id_filtered_luciphor_pep).join(ch_sdrf_config.luciphor_settings)

    output:
     set mzml_id, file("${id_file.baseName}_luciphor.idXML") into ptmt_in_id_luciphor
     file "*.log"

    when:
     params.enable_mod_localization

    script:
     def losses = params.luciphor_neutral_losses ? '-neutral_losses "${params.luciphor_neutral_losses}"' : ''
     def dec_mass = params.luciphor_decoy_mass ? '-decoy_mass "${params.luciphor_decoy_mass}"' : ''
     def dec_losses = params.luciphor_decoy_neutral_losses ? '-decoy_neutral_losses "${params.luciphor_decoy_neutral_losses}' : ''
     """
     LuciphorAdapter    -id ${id_file} \\
                        -in ${mzml_file} \\
                        -out ${id_file.baseName}_luciphor.idXML \\
                        -threads ${task.cpus} \\
                        -num_threads ${task.cpus} \\
                        -target_modifications ${params.mod_localization.tokenize(',').collect { "'${it}'" }.join(" ") } \\
                        -fragment_method ${frag_method} \\
                        ${losses} \\
                        ${dec_mass} \\
                        ${dec_losses} \\
                        -max_charge_state ${params.max_precursor_charge} \\
                        -max_peptide_length ${params.max_peptide_length} \\
                        -debug ${params.luciphor_debug} \\
                        > ${id_file.baseName}_luciphor.log
     """
                     //        -fragment_mass_tolerance ${} \\
                     //   -fragment_error_units ${} \\
}



process idmapper{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 tuple mzml_id, file(id_file_filter), file(consensusXML) from ptmt_in_id.mix(ptmt_in_id_luciphor).combine(id_files_consensusXML, by: 0)

	output:
	 file("${id_file_filter.baseName}_map.consensusXML") into id_map_to_merger

	script:
	 """
	 IDMapper -id ${id_file_filter} \\
	 		  -in ${consensusXML} \\
	 		  -threads ${task.cpus} \\
	 		  -rt_tolerance ${params.rt_tolerance} \\
	 		  -mz_tolerance ${params.mz_tolerance} \\
	 		  -mz_measure ${params.mz_measure} \\
	 		  -mz_reference ${params.mz_reference} \\
	 		  -debug 1 \\
	 		  -out ${id_file_filter.baseName}_map.consensusXML \\
	 		  > ${id_file_filter.baseName}_map.log
	 """
}


process file_merge{
	label 'process_medium'
	label 'process_single_thread'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 file(id_map) from id_map_to_merger.collect()

	output:
	 file("ID_mapper_merge.consensusXML") into id_merge_to_epi


	script:
	 """
	 FileMerger -in ${(id_map as List).join(' ')} \\
	 			-in_type consensusXML \\
	 			-annotate_file_origin  \\
	 			-append_method 'append_cols' \\
	 			-threads ${task.cpus} \\
	 			-debug 10 \\
	 			-out ID_mapper_merge.consensusXML \\
	 			> ID_mapper_merge.log
	 """
}


process epifany{
	label 'process_medium'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

	input:
	 file(consus_file) from id_merge_to_epi


	output:
	 file("${consus_file.baseName}_epi.consensusXML") into epi_idfilter

	 // expdes currently unused
	script:
	 """
	 Epifany -in ${consus_file} \\
	 		 -protein_fdr true \\
	 		 -threads ${task.cpus} \\
	 		 -debug 1 \\
       -algorithm:keep_best_PSM_only false \\
       -algorithm:update_PSM_probabilities false \\
	 		 -greedy_group_resolution ${params.greedy_group_resolution} \\
			 -algorithm:top_PSMs ${params.top_PSMs} \\
			 -out ${consus_file.baseName}_epi.consensusXML \\
			 > ${consus_file.baseName}_epi.log
	 """
}


process epi_filter{
  label 'process_very_low'
  label 'process_single_thread'

  publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

  input:
   file(consus_epi) from epi_idfilter


  output:
   file("${consus_epi.baseName}_filt.consensusXML") into conflict_res


  script:
  """
  IDFilter -in ${consus_epi} \\
           -out ${consus_epi.baseName}_filt.consensusXML \\
           -threads ${task.cpus} \\
           -score:prot ${params.protein_level_fdr_cutoff}
           > ${consus_epi.baseName}_idfilter.log
  """
}


process resolve_conflict{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
	publishDir "${params.outdir}/resolved_consensusXML", mode: 'copy', pattern: '*.consensusXML'

	input:
	 file(consus_epi_filt) from conflict_res


	output:
	 file "${consus_epi_filt.baseName}_resconf.consensusXML" into pro_quant


	script:
	"""
  IDConflictResolver -in ${consus_epi_filt} \\
						  -threads ${task.cpus} \\
						  -debug 1 \\
						  -resolve_between_features ${params.res_between_fet} \\
						  -out ${consus_epi_filt.baseName}_resconf.consensusXML \\
						  > ${consus_epi_filt.baseName}_resconf.log

	"""
}


process pro_quant{
	label 'process_medium'

	publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
	publishDir "${params.outdir}/proteomics_tmt", mode: 'copy'


	input:
	 file(epi_filt_resolve) from pro_quant
	 file pro_quant_exp from ch_pro_quant_exp


	output:
	 file "protein_out.csv" optional true into downstreams_tool_A
	 file "peptide_out.csv" into downstreams_tool_B
	 file "*.mzTab" optional true into out_mztab


	script:
    mztab = params.mztab_export ? "-mztab out.mzTab" : ""
    """
   	ProteinQuantifier -in ${epi_filt_resolve} \\
   				   	-design ${pro_quant_exp} \\
   				   	-out protein_out.csv \\
   				   	-peptide_out peptide_out.csv \\
             	${mztab} \\
   				   	-top ${params.top} \\
   				   	-average ${params.average} \\
   				   	-best_charge_and_fraction \\
   				   	-ratios \\
   				   	-threads ${task.cpus} \\
   				   	-consensus:normalize \\
   				   	-consensus:fix_peptides \\
   				   	-debug 100 \\
   				   	> pro_quant.log
   	"""
}


//TODO allow user config yml (as second arg to the script

process ptxqc {

    label 'process_low'
    label 'process_single_thread'

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'
    publishDir "${params.outdir}/ptxqc", mode: 'copy'

    when:
     params.enable_qc

    input:
     file mzTab from out_mztab

    output:
     file "*.html" into ch_ptxqc_report
     file "*.yaml"
     file "*.Rmd"
     file "*.pdf"
     file "*.txt"

    script:
     """
     ptxqc.R ${mzTab} > ptxqc.log
     """
}



//--------------------------------------------------------------- //
//---------------------- Nextflow specifics --------------------- //
//--------------------------------------------------------------- //


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-proteomicstmt-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/proteomicstmt Workflow Summary'
    section_href: 'https://github.com/nf-core/proteomicstmt'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    ThermoRawFileParser.sh --version &> v_thermorawfileparser.txt
    echo \$(FileConverter 2>&1) > v_fileconverter.txt || true
    echo \$(DecoyDatabase 2>&1) > v_decoydatabase.txt || true
    echo \$(CometAdapter 2>&1) > v_cometadapter.txt || true
    echo \$(comet 2>&1) > v_comet.txt || true
    echo \$(PeptideIndexer 2>&1) > v_peptideindexer.txt || true
    echo \$(PSMFeatureExtractor 2>&1) > v_psmfeatureextractor.txt || true
    echo \$(PercolatorAdapter 2>&1) > v_percolatoradapter.txt || true
    percolator -h &> v_percolator.txt
    echo \$(IDFilter 2>&1) > v_idfilter.txt || true
    echo \$(IDScoreSwitcher 2>&1) > v_idscoreswitcher.txt || true
    echo \$(FalseDiscoveryRate 2>&1) > v_falsediscoveryrate.txt || true
    echo \$(IDPosteriorErrorProbability 2>&1) > v_idposteriorerrorprobability.txt || true
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/proteomicstmt] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/proteomicstmt] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = ""
    try {
        if (workflow.success && ch_ptxqc_report.println()) {
            mqc_report = ch_ptxqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/proteomicstmt] Found multiple reports from process 'ptxqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
        else {
          mqc_report = ""
        }
    } catch (all) {
        log.warn "[nf-core/proteomicstmt] Could not attach PTXQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/proteomicstmt] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc != "" && mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/proteomicstmt] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/proteomicstmt]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/proteomicstmt]${c_red} Pipeline completed with errors${c_reset}-"
    }

}


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/proteomicstmt v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}


//--------------------------------------------------------------- //
//---------------------- Utility functions  --------------------- //
//--------------------------------------------------------------- //

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Check class of an Object for "List" type
boolean isCollectionOrArray(object) {
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}
