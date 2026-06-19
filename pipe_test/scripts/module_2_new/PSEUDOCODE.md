# structure_pipeline Pseudocode
Pseudocode for `module_2_new`, based on `README.md` and
`CODE_WALKTHROUGH.md`. The pipeline enriches LPMO protein inputs with UniProt
metadata, InterPro domains, LPMO feature parsing, optional BLAST
identification, and timestamped output files.
```text
PROGRAM structure_pipeline / module_2_new
MODULE STRUCTURE
    main_driver.py      = CLI, mode routing, enrichment orchestration, output
    input_handler.py    = FASTA, CAZy TSV, list, and characterized CSV parsing
    uniprot_client.py   = UniProt lookup, obsolete ID mapping, batch fetch, cache
    interpro_client.py  = InterPro domains, SignalP fallback, domain cleanup
    feature_parser.py   = signal peptide, LPMO core, CBMs, TM, H1 verification
    blast_client.py     = optional exact BLAST identification for unknown FASTA
MAIN
    Parse CLI:
        --mode = fasta | cazy | list | characterized
        --input, --output_dir, --cazy-family
        --allow-ncbi-fallback
        --allow-sequence-search
        --max-sequence-searches
    IF mode != cazy AND --cazy-family is missing:
        stop with CLI error
    timestamp = current time as YYYYMMDD_HHMMSS
    Create output folders:
        metadata/, sequences/, run/, cazy_raw/
    Define output files:
        metadata_expanded_TIMESTAMP.tsv
        all_sequences_TIMESTAMP.fasta
        failed_ids_TIMESTAMP.txt
        run_metadata_TIMESTAMP.json
    Create API clients with module-local JSON caches:
        UniProtClient(discovery_cache.json)
        InterProClient(interpro_cache.json)
    Initialize:
        successful_entries = []
        raw_fasta_entries = []
        failed_ids = []
        failed_id_set = set()
    SWITCH mode:
        cazy          -> run CAZy mode
        characterized -> run characterized mode
        fasta         -> run FASTA mode
        list          -> run list mode
    Write output files
    Write run metadata
MODE: CAZy
    IF input is CAZy family name or comma-separated family names:
        FOR each family:
            Download CAZy TSV
            Save raw TSV under cazy_raw/
            family_label = --cazy-family OR downloaded family
            Process downloaded file
    ELSE:
        family_label = --cazy-family OR family inferred from filename
        Process local CAZy TSV
    Process each CAZy TSV:
        CAZyHandler.process_cazy_file(file)
            Read tab-separated rows
            Split accessions into:
                ncbi_ids
                jgi_groups grouped by cleaned organism name
        UniProtClient.fetch_batch(ncbi_ids)
        FOR each returned NCBI-derived UniProt entry:
            enrich entry with family_label
        CAZyHandler.generate_jgi_queries(jgi_groups)
        FOR each JGI organism query:
            UniProtClient.search_by_query(query)
            FOR each returned JGI-matched UniProt entry:
                enrich entry with family_label
MODE: CHARACTERIZED
    CharacterizedHandler.process_characterized_file(input)
        Read semicolon CSV
        Require columns:
            Protein Name, EC#, Reference, Organism, GenBank, UniProt, PDB/3D
        Return:
            uni_ids
            genbank_ids
            rows_meta
    Split rows into direct UniProt rows, GenBank/RefSeq rows, and rows without IDs
    DEFINE choose_primary(id_list):
        Return first accession-like ID matching [A-Z0-9]{6,10}
        Otherwise return first ID
    DEFINE fetch_uniprot_ids(id_list):
        Fetch UniProt entries in batches
        FOR each entry:
            IF sequence is missing:
                mark ID as failed
            ELSE:
                fetch InterPro domains
                parse features
                store feature row by accession
                add accession to success set
                group identical sequences in seq_groups
    DEFINE map_non_uniprot_to_uniprot(id_token):
        IF ID looks like RefSeq:
            query UniProt with xref:RefSeq:ID
        ELSE:
            query UniProt with xref:EMBL-CDS:ID
        IF no result and ID has version suffix:
            retry without suffix
        Return first UniProt result or None
    DEFINE fetch_ncbi_fasta(genbank_acc):
        Fetch FASTA from NCBI efetch
        Parse sequence, organism, and protein name from header
        Return parsed values or empty values
    fetch_uniprot_ids(uni_ids)
    FOR each GenBank/RefSeq ID:
        mapped_entry = map_non_uniprot_to_uniprot(ID)
        IF mapped_entry exists:
            enrich mapped UniProt entry
            link original ID to mapped sequence group
        ELSE IF --allow-ncbi-fallback:
            fallback = fetch_ncbi_fasta(ID)
            IF fallback has sequence:
                create fallback metadata from NCBI + CSV row
                group fallback sequence
            ELSE:
                mark ID as failed
        ELSE:
            mark ID as failed
    FOR each identical sequence group:
        primary = choose_primary(group IDs)
        coaccessions = group IDs except primary
        create one FASTA record for the unique sequence
        add Sequence_Group and CoAccessions to metadata rows
    Detect CSV rows that produced multiple different sequences
    Store split-row count and row numbers for run_metadata.json
MODE: FASTA
    records = parse_fasta_file(input)
    FOR each FASTA record:
        parsed_id = parse_fasta_header(header)
        IF parsed_id is UniProt accession:
            add to known_ids
        ELSE:
            add header and sequence to unknown_headers
    UniProtClient.fetch_batch(known_ids)
    FOR each returned UniProt entry:
        enrich entry with CAZy_family = --cazy-family
    IF --allow-sequence-search:
        FOR each unknown sequence up to --max-sequence-searches:
            match = run_blast_search(sequence)
            IF exact 100 percent identity and coverage match exists:
                fetch matched UniProt entry
                enrich with Match_Status = BLAST_exact
            ELSE:
                mark unknown header as failed
        Log or mark unknown records skipped by the search limit
    ELSE:
        mark unknown headers as unresolved/failed
MODE: LIST
    Read input as plain text
    Strip whitespace and ignore empty lines
    Treat each remaining line as one UniProt accession
    UniProtClient.fetch_batch(accessions)
    FOR each returned UniProt entry:
        enrich entry with CAZy_family = --cazy-family
SHARED ENRICHMENT
    INPUT:
        UniProt entry or NCBI fallback metadata
        CAZy family label
        Match_Status
    Extract:
        accession, protein name, organism, EC number, sequence
    IF input is UniProt entry:
        InterProClient.fetch_domains(accession, sequence, signal_end)
            Use cache if available
            Query Pfam, CDD, SMART, PROSITE
            Parse domain coordinates
            Keep LPMO domains and binding modules
            Deduplicate overlaps by priority:
                Pfam > CDD > SMART > PROSITE
            Adjust LPMO start to first histidine after signal peptide
        feature_parser.parse_uniprot_features(entry, interpro_domains)
            Find signal peptide end
            Find transmembrane regions
            Find LPMO core start/end and provenance
            Find CBM binding modules
            Find H1 amino acid
            Set H1_Verified only when H1_AminoAcid is H
        IF UniProt signal peptide is missing:
            signal = InterProClient.fetch_signal_peptide(accession)
            IF signal is valid and inside sequence bounds:
                update Signal_End
                recompute H1_AminoAcid and H1_Verified
    ELSE IF input is NCBI fallback:
        Use available sequence, organism, protein name, and CSV metadata
        Fill unavailable UniProt/InterPro feature fields with defaults
    Append metadata row to successful_entries
    Append normalized FASTA record to raw_fasta_entries
UNIPROT CLIENT
    fetch_batch(ids):
        resolve obsolete/secondary IDs
        return cached entries when possible
        fetch missing IDs in batches
        recursively bisect failed batches to isolate bad IDs
        cache successful entries with atomic JSON write
    search_by_query(query):
        run custom UniProt search
        cache returned entries by primary accession
        return matching entries
INTERPRO CLIENT
    fetch_domains(accession, sequence, signal_end):
        return cached domains when possible
        fetch source-specific InterPro annotations
        deduplicate and adjust domains
        cache results with atomic JSON write
    fetch_signal_peptide(accession):
        return cached SignalP result when possible
        fetch SignalP annotation
        return corrected signal end or None
BLAST CLIENT
    run_blast_search(sequence):
        submit blastp search against SwissProt
        poll with timeout / max attempts
        accept exact 100 percent identity and coverage only
        return matched UniProt accession or None
OUTPUT
    Build pandas DataFrame from successful_entries
    Convert nullable integer columns where needed
    Write metadata TSV with:
        UniProt_ID, CAZy_family, InterPro_IDs, Protein_Name, Organism,
        EC_Number, Signal_End, Transmembrane_Regions, LPMO_Core_Start,
        LPMO_Core_End, Domain_Provenance, Binding_Modules, H1_Verified,
        H1_AminoAcid, Match_Status, Sequence_Group, CoAccessions
    Write FASTA:
        Standard modes:
            >UniProtID|ID|Organism|Protein_Name
        Characterized deduplicated mode:
            >UniProtIDs|ID1;ID2;...|Organism|Protein_Name
    Write failed_ids text file
    Write run_metadata JSON with:
        run_id, input, mode, counts, duration, output paths
        characterized-only sequence count and split-row details
ERROR HANDLING
    Catch KeyboardInterrupt separately
    Catch unexpected exceptions around main
    Log traceback
    Preserve partial outputs already written
    Exit with appropriate code
END PROGRAM
```
