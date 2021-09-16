rm Input_Data/ReferenceBased/* &> /dev/null
rm Input_Data/ReferenceFree/* &> /dev/null
rm -rf GeneratedFiles/out_spades_ &> /dev/null
rm GeneratedFiles/* &> /dev/null
rm src/SyntheticGenerator/Inputs/* &> /dev/null

# Clean Fasta Files Not compressed from DB
rm References/NCBI-Virus/GB_DB_VIRAL/*.fna &> /dev/null
rm References/NCBI-Virus/GB_DB_VIRAL_CDS/*.fna &> /dev/null
rm References/NCBI-Virus/GB_DB_VIRAL_RNA/*.fna &> /dev/null
rm References/NCBI-Archaea/GB_DB_ARCHAEA/*.fna &> /dev/null
rm References/NCBI-Archaea/GB_DB_ARCHAEA_CDS/*.fna &> /dev/null
rm References/NCBI-Archaea/GB_DB_ARCHAEA_RNA/*.fna &> /dev/null
rm References/NCBI-Bacteria/GB_DB_BACTERIA/*.fna &> /dev/null
rm References/NCBI-Bacteria/GB_DB_BACTERIA_CDS/*.fna &> /dev/null
rm References/NCBI-Bacteria/GB_DB_BACTERIA_RNA/*.fna &> /dev/null
rm References/NCBI-Fungi/GB_DB_FUNGI/*.fna &> /dev/null
rm References/NCBI-Fungi/GB_DB_FUNGI_CDS/*.fna &> /dev/null
rm References/NCBI-Fungi/GB_DB_FUNGI_RNA/*.fna &> /dev/null
rm References/NCBI-Mitochondrial/*.fna &> /dev/null
rm References/NCBI-Plant/GB_DB_PLANT/*.fna &> /dev/null
rm References/NCBI-Plant/GB_DB_PLANT_CDS/*.fna &> /dev/null
rm References/NCBI-Plant/GB_DB_PLANT_RNA/*.fna &> /dev/null
rm References/NCBI-Plastid/*.fna &> /dev/null
rm References/NCBI-Protozoa/GB_DB_PROTOZOA/*.fna &> /dev/null
rm References/NCBI-Protozoa/GB_DB_PROTOZOA_CDS/*.fna &> /dev/null
rm References/NCBI-Protozoa/GB_DB_PROTOZOA_RNA/*.fna &> /dev/null

# Clean Outputs
rm -rf Outputs/FalconNodes &> /dev/null
rm -rf Outputs/BlastnNodes &> /dev/null
rm Outputs/* &> /dev/null

# Clean Generated files from Input Sequences in Reference Free Analysis
rm -rf Input_Data/ReferenceFree/GeneratedORFs &> /dev/null
rm -rf Input_Data/ReferenceFree/GeneratedPredictorValues &> /dev/null

# Clean Result Files
rm Results/* &> /dev/null
rm -rf Results/falcon_seq &> /dev/null

# Clean cryfa files
rm Data_Security/Decrypted_Data/* &> /dev/null
rm Data_Security/Encrypt_Input/* &> /dev/null
rm Data_Security/Encrypted_Data/* &> /dev/null

# Clean Protein files
rm -rf ORFs/NodesProteins &> /dev/null

# Clean Reference Free Results
rm -rf Results/ReferenceFree &> /dev/null