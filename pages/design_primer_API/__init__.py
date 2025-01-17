# Copyright (c) 2023 Minniti Julien

# Portions of this software are based on TFinder, originally developed by Minniti Julien, under the MIT License.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files, to deal in the software
# without restriction, including without limitation the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies of the
# software, and to permit persons to whom the software is furnished to do so,
# subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import random
import re
import time
import xml.etree.ElementTree as ET

import primer3
import requests
from bs4 import BeautifulSoup
from tqdm import tqdm


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}


class NCBIdna:
    def __init__(self, gene_id, species=None, genome_version="Current", all_slice_forms=None):
        self.gene_id = gene_id
        self.species = species if species is not None else "human"
        self.all_slice_forms = True if all_slice_forms is True else False

    @staticmethod
    def XMNM_to_gene_ID(variant):
        global headers
        uids = f"https://www.ncbi.nlm.nih.gov/nuccore/{variant}"

        response = requests.get(uids, headers=headers)

        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')

            pattern = r"list_uids=(\d+)"
            matches = re.search(pattern, str(soup))

            if matches:
                entrez_id = matches.group(1)
            else:
                entrez_id = "UIDs not founds"
        else:
            entrez_id = "Error during process of retrieving UIDs"

        return entrez_id

    # Sequence extractor
    def find_sequences(self):
        time.sleep(1)
        if self.gene_id.startswith('XM_') or self.gene_id.startswith('NM_') or self.gene_id.startswith(
                'XR_') or self.gene_id.startswith('NR_'):
            entrez_id = NCBIdna.XMNM_to_gene_ID(self.gene_id)
            if entrez_id == 'UIDs not founds' or entrez_id == 'Error during process of retrieving UIDs':
                result_promoter = f'Please verify {self.gene_id} variant'
                return result_promoter
        else:
            if self.gene_id.isdigit():
                entrez_id = self.gene_id

            else:
                entrez_id, message = NCBIdna.convert_gene_to_entrez_id(self.gene_id, self.species)
                if entrez_id == "Error 200":
                    return entrez_id, message

        all_variants, message = NCBIdna.all_variant(entrez_id, self.all_slice_forms)
        if "Error 200" not in all_variants:
            return all_variants, message

    @staticmethod
    # Convert gene to ENTREZ_GENE_ID
    def convert_gene_to_entrez_id(gene_name, species):
        global headers

        while True:
            # Request for ENTREZ_GENE_ID
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="{gene_name}"[Gene%20Name]+AND+{species}[Organism]&retmode=json&rettype=xml'
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                response_data = response.json()

                if 'count' in response_data.get('esearchresult', {}):
                    if response_data['esearchresult']['count'] == '0':
                        print(bcolors.WARNING + f"Please verify if {gene_name} exist for {species}" + bcolors.ENDC)
                        return "Error 200", f"Please verify if {gene_name} exist for {species}"

                    else:
                        gene_id = response_data['esearchresult']['idlist'][0]
                        print(
                            bcolors.OKGREEN + f"Response 200: ID found for {species} {gene_name}: {gene_id}" + bcolors.ENDC)
                        return gene_id, bcolors.OKGREEN + f"Response 200: ID found for {species} {gene_name}: {gene_id}" + bcolors.ENDC
                else:
                    print(
                        bcolors.FAIL + f"Error 200: Issues for {species} {gene_name}, try again: {response.text}" + bcolors.ENDC)
                    time.sleep(random.uniform(0.25, 0.5))

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error 429: API rate limit exceeded during get ID of {species} {gene_name}, try again: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(bcolors.FAIL + f"Error {response.status_code}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get gene information
    def all_variant(entrez_id, all_slice_forms=False):
        global headers

        while True:
            url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
            response = requests.get(url2, headers=headers)

            if response.status_code == 200:
                response_data = response.json()
                try:
                    gene_info = response_data['result'][str(entrez_id)]
                    chraccver = gene_info['genomicinfo'][0]['chraccver']
                    print(
                        bcolors.OKGREEN + f"Response 200: Chromosome {chraccver} found for {entrez_id}: {response.text}" + bcolors.ENDC)
                    break

                except Exception as e:
                    print(
                        bcolors.WARNING + f"Response 200: Chromosome not found for {entrez_id}: {response.text} {e} {traceback.print_exc()}" + bcolors.ENDC)
                    all_variants = [("Error 200", None, None, None, None, None)]
                    print(
                        bcolors.WARNING + f"Response 200: Transcript not found(s) for {entrez_id}." + bcolors.ENDC)
                    return "Error 200", f"Transcript not found(s) for {entrez_id}."

            elif response.status_code == 429:
                print(
                    bcolors.ENDC + f"Error {response.status_code}: API rate limit exceeded during get chromosome of {entrez_id}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(
                    bcolors.ENDC + f"Error {response.status_code}: Error during get chromosome of {entrez_id}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={entrez_id}&retmode=xml"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                root = ET.fromstring(response.text)

                tv = []
                variants = []
                gene_name = []
                species_API = []
                chromosome = ""

                for elem in root.iter():
                    if elem.tag == "Gene-commentary_label":
                        if elem.text.startswith('transcript variant'):
                            if elem.text not in tv:
                                tv.append(elem.text)
                    if elem.tag == "Gene-commentary_accession":
                        if elem.text.startswith('NM_') or elem.text.startswith('XM_') or elem.text.startswith(
                                'NR_') or elem.text.startswith('XR_'):
                            if elem.text not in variants:
                                variants.append(elem.text)

                    elif elem.tag == "Org-ref_taxname":
                        species_API = elem.text

                    elif elem.tag == 'Gene-ref_locus':
                        gene_name = elem.text

                for elem in root.iter('Gene-commentary_accession'):
                    if elem.text.startswith('NC_'):
                        chromosome = elem.text
                        break

                def calc_exon(root, variants):
                    all_variants = {}
                    exon_coords = []
                    found_variant = False
                    k_found = False
                    orientation = ""

                    for variant in variants:
                        for elem in root.iter():
                            if elem.tag == "Gene-commentary_accession" and elem.text != variant:
                                if elem.text == chromosome:
                                    k_found = True
                                elif len(exon_coords) == 0:
                                    continue
                                else:
                                    break

                            if k_found and elem.tag == "Gene-commentary_accession":
                                found_variant = True if elem.text == variant else False
                            elif k_found and found_variant and elem.tag == "Seq-interval_from":
                                start = int(elem.text)
                            elif k_found and found_variant and elem.tag == "Seq-interval_to":
                                end = int(elem.text)
                                exon_coords.append((start, end))
                            elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                                orientation += elem.attrib.get("value")

                            elif elem.tag == "Org-ref_taxname":
                                species_API = elem.text

                            elif elem.tag == 'Gene-ref_locus':
                                gene_name = elem.text

                        if exon_coords:

                            if orientation == "minus":
                                exon_coords = [(end, start) for start, end in exon_coords]

                            first_exon_start = exon_coords[0][0]

                            normalized_exon_coords = [
                                (abs(start - first_exon_start), abs(end - first_exon_start)) for start, end in
                                exon_coords]

                            first_start = exon_coords[0][0]
                            last_end = exon_coords[-1][1]

                            sequence = NCBIdna.get_dna_sequence(entrez_id, chraccver, first_start, last_end)

                            all_variants[variant] = {
                                'entrez_id': entrez_id,
                                'gene_name': gene_name,
                                'chraccver': chraccver,
                                'exon_coords': exon_coords,
                                'normalized_exon_coords': normalized_exon_coords,
                                'species': species_API,
                                'sequence': sequence
                            }

                    if len(all_variants) > 0:
                        print(
                            bcolors.OKGREEN + f"Response 200: Transcript(s) found(s) for {entrez_id}: {all_variants}" + bcolors.ENDC)
                        return all_variants, f"Transcript(s) found(s) for {entrez_id}: {list(all_variants.keys())}"
                    else:
                        all_variants["Error 200"] = {
                            "entrez_id": f"Transcript not found for {entrez_id}.",
                            "gene_name": None,
                            "chraccver": None,
                            "exon_coords": None,
                            "normalized_exon_coords": None,
                            "species": None,
                            "sequence": None
                        }
                        print(
                            bcolors.WARNING + f"Error 200: Transcript not found(s) for {entrez_id}." + bcolors.ENDC)
                        return all_variants, f"Error 200: Transcript not found(s) for {entrez_id}."

                if all_slice_forms is True:
                    all_variants, message = calc_exon(root, variants)
                    return all_variants, message

                elif all_slice_forms is False:
                    if len(tv) > 0:
                        if "transcript variant 1" in tv:
                            associations = dict(zip(tv, variants))
                            variant = associations["transcript variant 1"]
                        else:
                            variant = variants[0]
                    elif len(tv) == 0 and len(variants) > 0:
                        variant = variants[0]
                    else:
                        variant = None

                    if variant is not None:
                        all_variants, message = calc_exon(root, [variant])
                        return all_variants, message

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error {response.status_code}: API rate limit exceeded while searching for {entrez_id} transcripts: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(
                    bcolors.FAIL + f"Error {response.status_code}: Error while searching for {entrez_id} transcripts: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get DNA sequence
    def get_dna_sequence(gene_name, chraccver, chrstart, chrstop):
        global headers

        if chrstop > chrstart:
            start = chrstart
            end = chrstop
        else:
            start = chrstop
            end = chrstart

        # Requête pour obtenir la séquence d'ADN
        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={chraccver}&from={start}&to={end}&rettype=fasta&retmode=text"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                dna_sequence = response.text.split('\n', 1)[1].replace('\n',
                                                                       '')

                if chrstop < chrstart:
                    sequence = NCBIdna.reverse_complement(dna_sequence)
                else:
                    sequence = dna_sequence

                print(f"Response 200: DNA sequence for {gene_name} extracted: {sequence}")
                return sequence

            elif response.status_code == 429:
                print(f"Error 429: API rate limit exceeded for DNA extraction of {gene_name}, try again.")
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(f"Error {response.status_code}: {response.text}")
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    def reverse_complement(dna_sequence):
        DNA_code = ["A", "T", "C", "G", "N", "a", "t", "c", "g", "n"]
        if not all(char in DNA_code for char in dna_sequence):
            isdna = 'Please use only A T G C'
            return isdna
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        reverse_sequence = dna_sequence[::-1].upper()
        complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
        return complement_sequence

    @staticmethod
    def design_primers(variant, gene_name, sequence, exons, nb_primers):
        try:
            simplified_sequence = "".join(sequence[start:end] for start, end in exons)

            primer3_input = {
                'SEQUENCE_ID': 'primer_in_exons',
                'SEQUENCE_TEMPLATE': simplified_sequence,
            }

            primer3_params = {
                'PRIMER_TASK': 'generic',  # Generic primer generation

                # Settings for primers (size and content)
                'PRIMER_OPT_SIZE': 20,  # Optimal primer size
                'PRIMER_MIN_SIZE': 16,  # Minimum primer size
                'PRIMER_MAX_SIZE': 24,  # Maximum primer size
                'PRIMER_OPT_TM': 60.0,  # Optimal melting temperature (°C)
                'PRIMER_MIN_TM': 53.0,  # Minimum melting temperature (°C)
                'PRIMER_MAX_TM': 63.0,  # Maximum melting temperature (°C)
                'PRIMER_MIN_GC': 45.0,  # Minimum GC percentage (%)
                'PRIMER_MAX_GC': 55.0,  # Maximum GC percentage (%)
                'PRIMER_GC_CLAMP': 0,  # GC clamping at end 3' (minimum number of G/C)
                'PRIMER_MAX_POLY_X': 5,  # Maximum number of repeated bases (ex: AAAAA)

                # Parameters for stability at the 3' end
                'PRIMER_MAX_END_STABILITY': 9.0,  # Maximum end stability 3'

                # Secondary alignment (Thermodynamic model)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 70.0,  # Bad template match (primer pairs)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TH_TMPL': 40.0,  # Bad match for single primer
                'PRIMER_MAX_SELF_ANY_TH': 45.0,  # Internal matching (all sites, thermodynamics)
                'PRIMER_MAX_SELF_END_TH': 35.0,  # Internal pairing (3' end, thermodynamic)
                'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,  # Primer pairing (all sites, thermodynamics)
                'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,  # Pairing between primers (3' end, thermodynamics)
                'PRIMER_MAX_HAIRPIN_TH': 24.0,  # Maximum free energy for hairpins

                # Secondary alignment (Old model)
                'PRIMER_MAX_TEMPLATE_MISPRIMING': 24.0,  # Bad template matching (primer pairs, classic)
                'PRIMER_MAX_TEMPLATE_MISPRIMING_TMPL': 12.0,  # Bad match for single primer (classic)
                'PRIMER_MAX_SELF_ANY': 8.0,  # Internal pairing (all sites, classic)
                'PRIMER_MAX_SELF_END': 3.0,  # Internal pairing (3' end, classic)
                'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,  # Primer pairing (all sites, classic)
                'PRIMER_PAIR_MAX_COMPL_END': 3.0,  # Pairing between primers (3' end, classic)

                # Search for secondary alignments
                'PRIMER_THERMODYNAMIC_ALIGNMENT': 1,  # Use thermodynamic model
                'PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT': 1 if len(sequence) < 10000 else 0,  # Also align with thermodynamic model (maybe slow)

                # General settings for pairs
                'PRIMER_NUM_RETURN': 1,  # Maximum number of pairs returned
                'PRIMER_PRODUCT_SIZE_RANGE': [[80, 250]],  # Product size range

                # Enable selection of internal hybridization oligos
                'PRIMER_PICK_INTERNAL_OLIGO': 0,  # 1 to enable, 0 to disable

                # Size parameters for internal oligos
                'PRIMER_INTERNAL_MIN_SIZE': 18,  # Minimum size
                'PRIMER_INTERNAL_OPT_SIZE': 20,  # Optimal size
                'PRIMER_INTERNAL_MAX_SIZE': 24,  # Maximum size

                # Melting temperature (Tm) for internal oligos
                'PRIMER_INTERNAL_MIN_TM': 57.0,  # Minimum Tm (°C)
                'PRIMER_INTERNAL_OPT_TM': 60.0,  # Optimal Tm (°C)
                'PRIMER_INTERNAL_MAX_TM': 63.0,  # Maximum Tm (°C)

                # GC percentage for internal oligos
                'PRIMER_INTERNAL_MIN_GC': 20.0,  # Minimum GC percentage
                'PRIMER_INTERNAL_OPT_GC_PERCENT': 50.0,  # Optimal GC percentage
                'PRIMER_INTERNAL_MAX_GC': 80.0,  # Maximum GC percentage

                # Concentration of monovalent cations (e.g.: Na+)
                'PRIMER_MONOVALENT_CATION_CONC': 50.0,  # In mM (default is 50 mM)

                # Concentration of divalent cations (e.g. Mg2+)
                'PRIMER_DIVALENT_CATION_CONC': 1.5,  # In mM (default is 1.5 mM)

                # Concentration of dNTPs
                'PRIMER_DNTP_CONC': 0.6,  # In mM (default is 0.8 mM)

                # Salt correction formula (using Santa Lucia 1998 formula)
                'PRIMER_SALT_CORRECTION': 1,  # 1 to use the Santa Lucia 1998 correction

                # Thermodynamic parameters (table of parameters to calculate Tm)
                'PRIMER_THERMODYNAMIC_PARAMETERS': 'SantaLucia1998',  # Use Santa Lucia 1998 table

                # Concentration of the oligonucleotide for the init
                'PRIMER_ANN_Oligo_CONC': 50.0,
            }

            primers = []

            exon_lengths = [end - start for start, end in exons]
            cumulative_lengths = [0] + list(NCBIdna.cumsum(exon_lengths))

            if len(exons) > 1:
                exon_pairs = [(i, j) for i in range(len(exons)) for j in range(i + 1, len(exons))]

            else:
                exon_pairs = [(0, 0)]

            seen_primers = set()

            with tqdm(total=nb_primers, desc=f"Generating primers for {variant} {gene_name}", unit="primer") as pbar:
                no_progress_count = 0
                max_no_progress = 2

                while len(primers) < nb_primers:
                    primer3_params['PRIMER_NUM_RETURN'] += 1
                    primers_found_in_iteration = False

                    for i, j in exon_pairs:
                        if len(primers) >= nb_primers:
                            break

                        exon1_start, exon1_end = exons[i]
                        exon2_start, exon2_end = exons[j]

                        simplified_start1 = cumulative_lengths[i]
                        simplified_end1 = simplified_start1 + (exon1_end - exon1_start)
                        simplified_start2 = cumulative_lengths[j]
                        simplified_end2 = simplified_start2 + (exon2_end - exon2_start)

                        product_size = simplified_end2 - simplified_start1

                        if 40 <= product_size:
                            primer3_input['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [
                                simplified_start1, simplified_end1 - simplified_start1,
                                simplified_start2, simplified_end2 - simplified_start2
                            ]

                            primer_results = primer3.bindings.design_primers(primer3_input, primer3_params)

                            # print(primer_results)

                            if 'PRIMER_PAIR_NUM_RETURNED' in primer_results and primer_results[
                                'PRIMER_PAIR_NUM_RETURNED'] > 0:
                                for k in range(primer_results['PRIMER_PAIR_NUM_RETURNED']):
                                    if len(primers) >= nb_primers:
                                        break

                                    left_key = f'PRIMER_LEFT_{k}_SEQUENCE'
                                    right_key = f'PRIMER_RIGHT_{k}_SEQUENCE'

                                    if left_key in primer_results and right_key in primer_results:
                                        left_seq = primer_results.get(left_key, 'N/A')
                                        right_seq = primer_results.get(right_key, 'N/A')

                                        primer_key = (left_seq, right_seq)
                                        if primer_key in seen_primers:
                                            continue

                                        left_position = primer_results.get(f'PRIMER_LEFT_{k}')[0]
                                        right_position = primer_results.get(f'PRIMER_RIGHT_{k}')[0]

                                        left_absolute = NCBIdna.convert_to_absolute(left_position, exons,
                                                                                    cumulative_lengths)
                                        right_absolute = NCBIdna.convert_to_absolute(right_position, exons,
                                                                                     cumulative_lengths)

                                        amplicon_size = right_position - left_position + 1
                                        amplicon_size_abs = right_absolute - left_absolute + 1

                                        primers.append({
                                            'left_primer': {
                                                'sequence': left_seq,
                                                'length': len(left_seq),
                                                'position': (left_position, left_position + len(left_seq)),
                                                'position_abs': (left_absolute, left_absolute + len(left_seq)),
                                                'tm': primer_results.get(f'PRIMER_LEFT_{k}_TM', 'N/A'),
                                                'gc_percent': primer_results.get(f'PRIMER_LEFT_{k}_GC_PERCENT', 'N/A'),
                                                'self_complementarity': primer_results.get(f'PRIMER_LEFT_{k}_SELF_ANY_TH',
                                                                                           'N/A'),
                                                'self_3prime_complementarity': primer_results.get(
                                                    f'PRIMER_LEFT_{k}_SELF_END_TH', 'N/A'),
                                                'exon_junction': None
                                            },
                                            'right_primer': {
                                                'sequence': right_seq,
                                                'length': len(right_seq),
                                                'position': (right_position - len(right_seq), right_position),
                                                'position_abs': (right_absolute - len(right_seq), right_absolute + len(right_seq)),
                                                'tm': primer_results.get(f'PRIMER_RIGHT_{k}_TM', 'N/A'),
                                                'gc_percent': primer_results.get(f'PRIMER_RIGHT_{k}_GC_PERCENT', 'N/A'),
                                                'self_complementarity': primer_results.get(f'PRIMER_RIGHT_{k}_SELF_ANY_TH',
                                                                                           'N/A'),
                                                'self_3prime_complementarity': primer_results.get(
                                                    f'PRIMER_RIGHT_{k}_SELF_END_TH', 'N/A'),
                                                'template_strand': 'Minus',
                                                'exon_junction': None
                                            },
                                            'amplicon_size': amplicon_size,
                                            'amplicon_size_abs': amplicon_size_abs
                                        })

                                        pbar.update(1)
                                        seen_primers.add(primer_key)
                                        primers_found_in_iteration = True

                    if primers_found_in_iteration is False:
                        no_progress_count += 1
                    else:
                        no_progress_count = 0

                    if no_progress_count >= max_no_progress:
                        print("Breaking the loop: No new primers found after 5 iterations.")
                        break

            return primers

        except Exception as e:
            print(e)


    @staticmethod
    def cumsum(iterable):
        total = 0
        for value in iterable:
            total += value
            yield total

    @staticmethod
    def convert_to_absolute(position, exons, cumulative_lengths):
        for i, (start, end) in enumerate(exons):
            if position < cumulative_lengths[i + 1]:
                offset = position - cumulative_lengths[i]
                return start + offset
        return -1
