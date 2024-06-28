import argparse
import requests
import json


def download_germline_set(species, locus, version='latest', germline_set_name=None, file_format='AIRRC-JSON', url='https://ogrdb.airr-community.org/api_v2', prefix=None):
    """
    Downloads the specified germline set using the OGRDB API.

    Args:
        species (str): The species name (e.g., "Human").
        locus (str): The locus (e.g., IGH, IGK, IGL, TRA, TRB, TRD).
        version (str): The specific version to download, default is 'latest'.
        germline_set_name (str): The name of the germline set, default is None.
        file_format (str): The format to download (AIRRC-JSON, SINGLE-FG, SINGLE-FU).
        url (str): The base URL for the API, default is 'https://ogrdb.airr-community.org/api_v1'.
        prefix (str): The prefix for filenames, default is None.
    """
    found_species = find_species(url, species)
    if found_species:
        germline_set_id = find_available_sets(url, species, locus, germline_set_name)
        if germline_set_id is not None:
            data = request_set(url, germline_set_id, version, species, germline_set_name, locus, file_format)
            if not data:
                return
            if prefix is None:
                prefix = species.replace(' ', '_')
            if file_format == 'AIRRC-JSON':
                output_file = f"{prefix}.json"
                write_json_to_disk(data, output_file)
            else:
                if file_format == 'SINGLE-FG':
                    output_file = f"{prefix}_gapped.fasta"
                elif file_format == 'SINGLE-FU':
                    output_file = f"{prefix}_gapped.fasta"

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
    versions_list = get_set_versions(url, germline_set_id)
    if version == "latest" or version in versions_list:
        if format == 'AIRRC-JSON':
            api_url = f"{url}/germline/set/{germline_set_id}/{version}"
            response = requests.get(api_url)
            response.raise_for_status()
            return response.json()
        else:
            if format == "SINGLE-FG":
                api_url = f"{url}/germline/set/{germline_set_id}/{version}/gapped"
            elif format == "SINGLE-FU":
                api_url = f"{url}/germline/set/{germline_set_id}/{version}/ungapped"

            response = requests.get(api_url)
            response.raise_for_status()
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
    response = requests.get(api_url)
    response.raise_for_status()
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
    response = requests.get(api_url)
    response.raise_for_status()
    sets_list = response.json().get('germline_species')
    found = False
    found_count = 0
    available_sets = []
    germline_set_id = None
    for set in sets_list:
        if set.get('locus') == locus:
            found = True
            found_count += 1
            available_sets.append(set.get('germline_set_name'))
            germline_set_id = set.get('germline_set_id')
    
    if not found:
        available_locus = []
        for obj in sets_list:
            if obj.get('locus') not in available_locus:
                available_locus.append(obj.get('locus'))  

        print(f'''No sets for locus {locus} in species {species_id}. Available loci:
{', '.join(available_locus)}''')

    elif found_count > 1 and germline_set_name is None:
        germline_set_id = None
        print(f'''Multiple sets are avaliable for species {species_id} locus {locus}. 
Please specify the set name. Available set names:
{', '.join(available_sets)}''')

    else:
        print(f'{germline_set_id}')

    return germline_set_id


def find_species(url, species):
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
    response = requests.get(api_url)
    response.raise_for_status()
    species_list = response.json().get('species')
    found = False
    for object in species_list:
        if object.get('label') == species:
            found = True

    if not found:
        available_species = ', '.join(obj.get('label') for obj in species_list)
        print(f'''Error: species not found. Available species:
{available_species}''')
  
    else:
        print(f'{species}')

    return found


def get_parser():
    parser = argparse.ArgumentParser(description='Download germline sets from the Open Germline Receptor Database (OGRDB)')
    parser.add_argument('species', type=str, help='Species (e.g. "Human")')
    parser.add_argument('locus', type=str, help='Locus (IGH, IGK, IGL, TRA, TRB, TRD, TRD)')
    parser.add_argument('-n', '--name', type=str, help='germline set name (the utility will attempt to determine the name, if none is specified)')
    parser.add_argument('-v', '--version', type=str, default='latest', help='Specific version to download, otherwise the latest version will be downloaded')
    parser.add_argument('-f', '--format', type=str, default='AIRRC-JSON', choices=['AIRRC-JSON', 'SINGLE-FG', 'SINGLE-FU'], help='Format to download')
    parser.add_argument('-u', '--url', type=str, default='https://ogrdb.airr-community.org/api_v1', help='URL to use')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix for filenames. Default prefix is species (with _ substituted for space). If PREFIX is NONE, no prefix will be used for multi files, and the default prefix will be used for single files.')
    return parser


def main():
    args = get_parser().parse_args()
    download_germline_set(args.species, args.locus, args.version, args.name, args.format, args.url, args.prefix)


if __name__ == "__main__":
    main()