import argparse
import requests
import json
from Bio import SeqIO
from io import StringIO
from receptor_utils import simple_bio_seq as simple
from receptor_utils import aux_formats


def do_request(url):
    """
    Makes a GET request to the specified URL.

    Args:
        url (str): The URL to request.

    Returns:
        requests.Response: The response object.
    """

    # handle errors
    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.HTTPError as err:
        print(f"HTTP error occurred: {err} {response.text}")
        quit()
    except requests.exceptions.RequestException as err:
        print(f"Request error occurred: {err}")
        quit()

    return response


def download_germline_set(species, locus, version='latest', germline_set_name=None, file_format='AIRRC-JSON', url='https://ogrdb.airr-community.org/api_v2', prefix=None):
    """
    Downloads the specified germline set using the OGRDB API.

    Args:
        species (str): The species name (e.g., "Human").
        locus (str): The locus (e.g., IGH, IGK, IGL, TRA, TRB, TRD).
        version (str): The specific version to download, default is 'latest'.
        germline_set_name (str): The name of the germline set, default is None.
        file_format (str): The format to download (AIRRC-JSON, SINGLE-FG, SINGLE-FU).
        url (str): The base URL for the API, default is 'https://ogrdb.airr-community.org/api_v2'.
        prefix (str): The prefix for filenames, default is None.
    """
    species_id = find_species_id(url, species)
    if species_id:
        germline_set_id = find_available_sets(url, species_id, locus, germline_set_name)
        if germline_set_id is not None:
            if prefix is None:
                if germline_set_name:
                    prefix = species.replace(' ', '_') + '_' + germline_set_name.replace(' ', '_').replace('/', '_') + '_' + locus + '_'
                else:
                    prefix = species.replace(' ', '_') + '_' + locus + '_'
            elif prefix == 'NONE':
                prefix = ''

            if 'MULTI' in file_format:
                data_u = request_set(url, germline_set_id, version, species, germline_set_name, locus, 'SINGLE-FU')
                data_g = request_set(url, germline_set_id, version, species, germline_set_name, locus, 'SINGLE-FG')

                if not data_u or not data_g:
                    return
                
                data_u = SeqIO.parse(StringIO(data_u.text), 'fasta')
                data_u = {rec.id: str(rec.seq).upper() for rec in data_u}
                data_g = SeqIO.parse(StringIO(data_g.text), 'fasta')
                data_g = {rec.id: str(rec.seq).upper() for rec in data_g}
                files_written = []
                for chain in ['V', 'D', 'J', 'C']:
                    if chain in ['V', 'J']:
                        seqs = {k: v for k, v in data_u.items() if k[3] == chain}
                    elif chain == 'D':
                        seqs = {k: v for k, v in data_u.items() if k[3] == 'D' and len(v) < 100}    # Hopefully this gives good enough differentiation
                    else:
                        seqs = {k: v for k, v in data_g.items() if k[3] not in ['V', 'D', 'J'] or (k[3] == 'D' and len(v) >= 100)}
                    if seqs: 
                        simple.write_fasta(f"{prefix}{chain}.fasta", seqs)
                        files_written.append(f"{prefix}{chain}.fasta")
                seqs = {k: v for k, v in data_g.items() if k[3] == 'V'}
                if seqs:
                    simple.write_fasta(f"{prefix}V_gapped.fasta", seqs)
                    files_written.append(f"{prefix}V_gapped.fasta")
                print(f"FASTA files saved to {', '.join(files_written)}")

                if file_format == 'MULTI-IGBLAST':
                    data = request_set(url, germline_set_id, version, species, germline_set_name, locus, 'AIRRC-JSON')
                    write_igblast_to_disk(data, locus, prefix)
            else:
                data = request_set(url, germline_set_id, version, species, germline_set_name, locus, file_format)
                if not data:
                    return
                if file_format == 'AIRRC-JSON':
                    if prefix[-1] == '_':
                        prefix = prefix[:-1]
                    output_file = f"{prefix}.json"
                    write_json_to_disk(data, output_file)
                elif file_format == 'SINGLE-FG':
                    output_file = f"{prefix}gapped.fasta"
                    write_fasta_to_disk(data.text, output_file)
                elif file_format == 'SINGLE-FU':
                    output_file = f"{prefix}ungapped.fasta"
                    write_fasta_to_disk(data.text, output_file)



def write_fasta_to_disk(data, output_file):
    """
    Writes FASTA data to a file.

    Args:
        data (str): The FASTA data to write.
        output_file (str): The name of the output file.
    """
    with open(output_file, 'w') as f:
        f.write(data)
    
    print(f"FASTA file saved to {output_file}")


def write_json_to_disk(data, output_file):
    """
    Writes JSON data to a file.

    Args:
        data (dict): The JSON data to write.
        output_file (str): The name of the output file.
    """
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"JSON file saved to {output_file}")


def write_igblast_to_disk(data, locus, prefix):
    """
    Writes IGBLAST data to a file.

    Args:
        data (str): The JSON data structure.
        output_file (str): The name of the output file.
    """
    p = prefix
    if p[-1] == '_':
        p = p[:-1]
    if aux_formats.ndm_from_json(data, f"{p}.ndm", f"V{locus[2]}", 'IMGT'):
        print(f"IGBLAST ndm file saved to {p}.ndm")
    else:
        print("No V genes found in reference set. Skipping ndm file creation.")

    if aux_formats.aux_from_json(data, f"{p}.aux", False):
        print(f"IGBLAST aux file saved to {p}.aux")
    else:
        print("No J genes found in reference set. Skipping aux file creation.")


def request_set(url, germline_set_id, version, species, germline_set_name, locus, format):
    """
    Requests the germline set data from the API.

    Args:
        url (str): The base URL for the API.
        germline_set_id (str): The ID of the germline set.
        version (str): The version to download.
        species (str): The species name.
        germline_set_name (str): The name of the germline set.
        file_format (str): The format to download (AIRRC-JSON, SINGLE-FG, SINGLE-FU).

    Returns:
        dict or requests.Response: The requested data in JSON or FASTA format.
    """

    if version != 'latest':
        versions_list = get_set_versions(url, germline_set_id)
    else:
        versions_list = []

    ex = '_ex' if species == 'Homo sapiens' else ''

    if version == "latest" or version in versions_list:
        if format == 'AIRRC-JSON':
            api_url = f"{url}/germline/set/{germline_set_id}/{version}"
            response = do_request(api_url)
            return response.json()
        else:
            if format == "SINGLE-FG":
                api_url = f"{url}/germline/set/{germline_set_id}/{version}/gapped{ex}"
            elif format == "SINGLE-FU":
                api_url = f"{url}/germline/set/{germline_set_id}/{version}/ungapped{ex}"

            response = do_request(api_url)
            return response 
    else:
        gs = germline_set_name if germline_set_name else ""
        loc = locus if locus else ""
        print(f"""Version {version} of set {species} {gs} {loc} is not available. Available versions:
{', '.join(versions_list)}""")
        return None


def get_set_versions(url, germline_set_id):
    """
    Gets the available versions of a germline set.

    Args:
        url (str): The base URL for the API.
        germline_set_id (str): The ID of the germline set.

    Returns:
        list: A list of available versions.
    """
    api_url = f"{url}/germline/set/{germline_set_id}/versions"
    response = do_request(api_url)
    versions_list = [str(x) for x in response.json()['versions']]
    return versions_list


def find_available_sets(url, species_id, locus, germline_set_name):
    """
    Finds available germline sets for a species and locus.

    Args:
        url (str): The base URL for the API.
        species_id (str): The species ID.
        locus (str): The locus (e.g., IGH, IGK, IGL, TRA, TRB, TRD).
        germline_set_name (str): The name of the germline set.

    Returns:
        str: The ID of the found germline set, or None if not found.
    """
    api_url = f"{url}/germline/sets/{species_id}"
    response = do_request(api_url)
    sets_list = response.json().get('germline_species')
    found = False
    locus_found = False
    found_count = 0
    available_sets = []
    available_versions = []
    germline_set_id = None
    for curr_set in sets_list:
        if curr_set.get('locus') == locus:
            locus_found = True
            available_sets.append(curr_set['germline_set_name'])
        if curr_set.get('locus') == locus and (germline_set_name is None or germline_set_name == curr_set['germline_set_name']):
            found = True
            found_count += 1
            available_versions.append(curr_set['germline_set_id'])
            germline_set_id = curr_set['germline_set_id']
    
    if not locus_found:
        available_locus = []
        for obj in sets_list:
            if obj.get('locus') not in available_locus:
                available_locus.append(obj.get('locus'))  

        print(f'''No sets for locus {locus} in species {species_id}. Available loci:
{', '.join(available_locus)}''')
        
    if not found and germline_set_name:
        print(f'''Error: set {germline_set_name} not found for species {species_id} locus {locus}. Available sets:
{', '.join(available_sets)}''')
        
    elif found_count > 1 and germline_set_name is None:
        if len(set(available_sets)) >= 1:
            germline_set_id = None
            print(f'''Multiple sets are avaliable for species {species_id} locus {locus}. 
Please specify the set name. Available set names:
{', '.join(available_sets)}''')
        else:
            return 'latest'
    else:
        print(f'{germline_set_id}')

    return germline_set_id


def find_species_id(url, species):
    """
    Finds if a species is available in the API.

    Args:
        url (str): The base URL for the API.
        species (str): The species name.

    Returns:
        bool: True if the species is found, False otherwise.
    """
    api_url = f"{url}/germline/species"
    print(api_url)
    response = do_request(api_url)

    species_list = response.json().get('species')
    id = None
    for object in species_list:
        if object.get('label') == species:
            id = object.get('id')

    if not id:
        available_species = ', '.join(obj.get('label') for obj in species_list)
        print(f'''Error: species not found. Available species:
{available_species}''')

    else:
        print(f'{species}: {id}')

    return id


def get_parser():
    parser = argparse.ArgumentParser(description='Download germline sets from the Open Germline Receptor Database (OGRDB)')
    parser.add_argument('species', type=str, help='Species (e.g. "Homo sapiens")')
    parser.add_argument('locus', type=str, help='Locus (IGH, IGK, IGL, TRA, TRB, TRD, TRD)')
    parser.add_argument('-n', '--name', type=str, help='germline set name (the utility will attempt to determine the name, if none is specified)')
    parser.add_argument('-v', '--version', type=str, default='latest', help='Specific version to download, otherwise the latest version will be downloaded')
    parser.add_argument('-f', '--format', type=str, default='AIRRC-JSON', choices=['AIRRC-JSON', 'SINGLE-FG', 'SINGLE-FU', 'MULTI-F', 'MULTI-IGBLAST'], help='Format to download')
    parser.add_argument('-u', '--url', type=str, default='https://ogrdb.airr-community.org/api_v2', help='URL to use')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix for filenames. Default prefix is species_locus (with _ substituted for space). If PREFIX is NONE, no prefix will be used for multi files, and the default prefix will be used for single files.')
    return parser


def main():
    args = get_parser().parse_args()
    download_germline_set(args.species, args.locus, args.version, args.name, args.format, args.url, args.prefix)


if __name__ == "__main__":
    main()
