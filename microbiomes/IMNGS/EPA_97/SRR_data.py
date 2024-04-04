from Bio import Entrez

# Replace 'your_email@example.com' with your email address
Entrez.email = 'abonacolta@gmail.com'

def fetch_metadata_for_srr(srr_number):
    handle = Entrez.esearch(db='sra', term=srr_number)
    record = Entrez.read(handle)
    if record['IdList']:
        sra_id = record['IdList'][0]
        handle = Entrez.efetch(db='sra', id=sra_id, rettype='runinfo', retmode='text')
        metadata = handle.read().decode('utf-8')  # Decode the bytes to a string
        return metadata

def fetch_metadata_for_srr_list(file_path):
    with open(file_path, 'r') as file:
        srr_list = [line.strip() for line in file if line.strip()]

    all_metadata = []

    for srr_number in srr_list:
        metadata = fetch_metadata_for_srr(srr_number)
        all_metadata.append(metadata)

    return all_metadata

def save_metadata_to_file(metadata_list, output_file):
    with open(output_file, 'w') as file:
        for metadata in metadata_list:
            file.write(metadata)

# Replace 'SRR_list.txt' with the path to your text file containing SRR numbers
# Replace 'all_metadata.txt' with the desired output file
srr_list_file = 'SRRs.txt'
output_metadata_file = 'all_metadata.txt'

metadata_list = fetch_metadata_for_srr_list(srr_list_file)
save_metadata_to_file(metadata_list, output_metadata_file)
