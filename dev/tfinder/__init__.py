# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of TFinder and associated documentation files, to deal
# in TFinder without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of TFinder, and to permit persons to whom TFinder is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of TFinder.

# TFINDER IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH TFINDER OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import random
import re
import time
import xml.etree.ElementTree as ET

import altair as alt
import logomaker
import numpy as np
import pandas as pd
import requests
import streamlit as st
from bs4 import BeautifulSoup
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import StandardScaler
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
        self.genome_version = genome_version if genome_version is not None else "Current"
        self.all_slice_forms = True if all_slice_forms is True else False

    @staticmethod
    def XMNM_to_gene_ID(variant):
        uids = f"https://www.ncbi.nlm.nih.gov/nuccore/{variant}"

        response = requests.get(uids)

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

    @staticmethod
    # Analyse if gene is available
    def analyse_gene(gene_id):
        disponibility_list = ['ID', 'Human', 'Mouse', 'Rat', 'Drosophila', 'Zebrafish']
        time.sleep(0.25)
        gene_analyse = [gene_id]
        for species_test in disponibility_list:
            if not gene_id.isdigit():
                if species_test == 'ID':
                    gene_analyse.append('n.d')
                else:
                    time.sleep(0.5)
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={gene_id}[Gene%20Name]+AND+{species_test}[Organism]&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if response_data['esearchresult']['count'] != '0':
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

            if gene_id.isdigit():
                if species_test != 'ID':
                    gene_analyse.append('n.d')
                else:
                    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={gene_id}&retmode=json&rettype=xml"
                    response = requests.get(url)

                    if response.status_code == 200:
                        response_data = response.json()

                        if 'chraccver' in str(response_data):
                            gene_analyse.append("✅")
                        else:
                            gene_analyse.append("❌")

        return gene_analyse

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
                variant, gene_name, title, chraccver, chrstart, chrstop, strand, species_API = NCBIdna.get_variant_info(
                    entrez_id,
                    self.gene_id)
        else:
            if self.gene_id.isdigit():
                entrez_id = self.gene_id

            else:
                entrez_id, message = NCBIdna.convert_gene_to_entrez_id(self.gene_id, self.species)
                if entrez_id == "Error 200":
                    return entrez_id, message

            variant, gene_name, title, chraccver, chrstart, chrstop, strand, species_API, message = NCBIdna.get_gene_info(
                entrez_id, self.genome_version, gene_name_error=self.gene_id)
            if variant == "Error 200":
                return variant, message

            if self.all_slice_forms:
                all_variants, message = NCBIdna.all_variant(entrez_id)
                if all_variants == "Error 200":
                    all_variants = [(variant, gene_name, chraccver, exon_coords, normalized_exon_coords, species_API)]

        prom_term = self.prom_term.lower()
        if prom_term not in ['promoter', 'terminator']:
            result_promoter = f"'{self.prom_term}' not valid. Please use 'Promoter' or 'Terminator'."
            return result_promoter, "OK"

        if isinstance(self.upstream, int) and isinstance(self.downstream, int):
            upstream = int(self.upstream)
            downstream = int(self.downstream)
        else:
            result_window = f'Upstream {self.upstream} and Downstream {self.downstream} must be integer'
            return result_window, "OK"

        if not self.all_slice_forms or self.all_slice_forms and self.gene_id.startswith(
                'XM_') or self.gene_id.startswith('NM_') or self.gene_id.startswith('XR_') or self.gene_id.startswith(
            'NR_'):
            dna_sequence = NCBIdna.get_dna_sequence(gene_name, chraccver, chrstart, chrstop)

            dna_sequence = f">{variant} {gene_name} | {title} {chraccver} | Strand: {strand} | TSS (on chromosome): {chrstart + 1}\n"

            return dna_sequence, "OK"

        elif self.all_slice_forms:
            result_compil = []
            for variant, gene_name, chraccver, exon_coords, _, species_API in all_variants:
                chrstart = exon_coords[0][0]
                chrstop = exon_coords[-1][1]
                dna_sequence = NCBIdna.get_dna_sequence(gene_name, chraccver, chrstart, chrstop)
                results = f">{variant} {gene_name} | {title} | {chraccver} | Strand: {strand} | TSS (on chromosome): {chrstart + 1}\n"

                result_compil.append(results)

            result_output = "\n".join(result_compil)
            return result_output, "OK"

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
    def get_gene_info(entrez_id, genome_version="Current", from_id=True, gene_name_error=None):
        global headers

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
            response = requests.get(url, headers=headers)

            if response.status_code == 200:
                response_data = response.json()
                try:
                    gene_info = response_data['result'][str(entrez_id)]
                    gene_name = gene_info['name']
                    title, chraccver, chrstart, chrstop = NCBIdna.extract_genomic_info(entrez_id, response_data,
                                                                                       genome_version)
                    species_API = gene_info['organism']['scientificname']

                    if from_id:
                        variant = NCBIdna.all_variant(entrez_id, from_id=True)
                        if not variant:
                            variant = gene_name

                    strand = "plus" if chrstart < chrstop else "minus"

                    print(bcolors.OKGREEN + f"Response 200: Info for {entrez_id} {variant} {gene_name} retrieved."
                                            f"Entrez_ID: {entrez_id} | Gene name: {gene_name} | Genome Assembly: {title} | ChrAccVer: {chraccver}"
                                            f"ChrStart/ChrStop: {chrstart}/{chrstop} | Strand: {strand} | Species: {species_API}" + bcolors.ENDC)

                    return variant if variant else None, gene_name, title, chraccver, chrstart, chrstop, strand, species_API, (
                        f"Response 200: Info for {entrez_id} {variant} {gene_name} retrieved."
                        f"Entrez_ID: {entrez_id} | Gene name: {gene_name} | Genome Assembly: {title} | ChrAccVer: {chraccver}"
                        f"ChrStart/ChrStop: {chrstart}/{chrstop} | Strand: {strand} | Species: {species_API}")
                except Exception as e:
                    print(
                        bcolors.WARNING + f"Response 200: Info for {gene_name_error} {entrez_id} not found: {e} {traceback.print_exc()}" + bcolors.ENDC)
                    return "Error 200", None, None, None, None, None, None, None, f"Info for {gene_name_error} {entrez_id} not found."

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"Error 429: API rate limit exceeded during get {entrez_id} info, try again: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(bcolors.FAIL + f"Error {response.status_code}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    def extract_genomic_info(gene_id, gene_info, genome_version):
        if gene_info and 'result' in gene_info and gene_id in gene_info['result']:
            accession_dict = {}
            gene_details = gene_info['result'][gene_id]

            time.sleep(1)

            location_hist = gene_details.get('locationhist', [])
            if len(location_hist) == 0:
                location_hist = gene_details.get('genomicinfo', [])
            for loc in location_hist:
                nc_accver = loc.get('chraccver')
                chrstart = loc.get('chrstart')
                chrstop = loc.get('chrstop')

                if nc_accver:
                    base_accession = nc_accver
                    if base_accession not in accession_dict:
                        accession_dict[base_accession] = (chrstart, chrstop)
                    else:
                        existing_start, existing_stop = accession_dict[base_accession]
                        accession_dict[base_accession] = (min(existing_start, chrstart), max(existing_stop, chrstop))

            nc_dict = accession_dict

            nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                       base_accver.startswith(("NC_", "NT_"))}

            if nc_dict:
                first_base = next(iter(nc_dict)).split('.')[0]

                nc_dict = {base_accver: (chrstart, chrstop) for base_accver, (chrstart, chrstop) in nc_dict.items() if
                           base_accver.split('.')[0] == first_base}

            max_version = -1
            max_accver = None
            max_coords = None
            min_version = float('inf')
            min_accver = None
            min_coords = None

            for base_accver in nc_dict.keys():
                version = int(base_accver.split('.')[1])

                if version > max_version:
                    max_version = version
                    max_accver = base_accver
                    max_coords = nc_dict[base_accver]

                if version < min_version:
                    min_version = version
                    min_accver = base_accver
                    min_coords = nc_dict[base_accver]

            if genome_version != "Current":
                title = NCBIdna.fetch_nc_info(min_accver)
                return title, min_accver, min_coords[0], min_coords[1]
            else:
                title = NCBIdna.fetch_nc_info(max_accver)
                return title, max_accver, max_coords[0], max_coords[1]

    @staticmethod
    def fetch_nc_info(nc_accver):
        global headers

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={nc_accver}&retmode=json"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                nc_info = response.json()
                try:
                    uid = nc_info['result']['uids'][0]
                    title = nc_info['result'][uid]['title']
                    return title
                except Exception as e:
                    time.sleep(random.uniform(0.25, 0.5))
            else:
                time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get gene information
    def get_variant_info(entrez_id, variant):
        global headers
        variant = variant.split(".")[0]

        while True:
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={entrez_id}&retmode=xml"
            response = requests.get(url, headers=headers)
            if response.status_code == 200:
                root = ET.fromstring(response.text)

                chromosome = ""
                found_variant = False
                k_found = False
                start_coords = []
                end_coords = []
                orientation = ""
                break

            elif response.status_code == 429:
                print(
                    bcolors.FAIL + f"1 Error 429: API rate limit exceeded during get {entrez_id} {variant} info, try again: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))
            else:
                print(bcolors.FAIL + f"Error {response.status_code}: {response.text}" + bcolors.ENDC)
                time.sleep(random.uniform(0.25, 0.5))

        if response.status_code == 200:
            for elem in root.iter('Gene-commentary_accession'):
                if elem.text.startswith('NC_'):
                    chromosome = elem.text
                    break

            for elem in root.iter():
                if elem.tag == "Gene-commentary_accession" and elem.text != variant:
                    if elem.text == chromosome:
                        k_found = True
                    elif len(start_coords) < 2 and len(end_coords) < 2:
                        continue
                    else:
                        break

                if k_found and elem.tag == "Gene-commentary_accession":
                    found_variant = True if elem.text == variant else False
                elif k_found and found_variant and elem.tag == "Seq-interval_from":
                    start_coords.append(elem.text)
                elif k_found and found_variant and elem.tag == "Seq-interval_to":
                    end_coords.append(elem.text)
                elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                    orientation += elem.attrib.get("value")

                elif elem.tag == "Org-ref_taxname":
                    species_API = elem.text

                elif elem.tag == 'Gene-ref_locus':
                    gene_name = elem.text

            if orientation != "minus":
                chrstart = int(start_coords[0]) + 1
                chrstop = int(end_coords[-1]) + 1
            else:
                chrstart = int(end_coords[0]) + 1
                chrstop = int(start_coords[-1]) + 1

            strand = "plus" if chrstart < chrstop else "minus"

            while True:
                url2 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id={entrez_id}&retmode=json&rettype=xml"
                response = requests.get(url2, headers=headers)
                if response.status_code == 200:
                    response_data = response.json()
                    if 'result' in response_data and str(entrez_id) in response_data['result']:
                        gene_info = response_data['result'][str(entrez_id)]
                        if 'chraccver' in str(gene_info):
                            chraccver = gene_info['genomicinfo'][0]['chraccver']
                            title = NCBIdna.fetch_nc_info(chraccver)
                    else:
                        gene_name = 'Bad ID'

                    print(bcolors.OKGREEN + f"Response 200: Info for {entrez_id} {variant} {gene_name} retrieved. "
                                            f"Entrez_ID: {entrez_id} | Gene name: {gene_name} | ChrAccVer: {chraccver}"
                                            f"ChrStart/ChrStop: {chrstart}/{chrstop} | Species: {species_API}" + bcolors.ENDC)
                    return variant, gene_name, title, chraccver, chrstart, chrstop, strand, species_API

                elif response.status_code == 429:
                    print(
                        bcolors.FAIL + f"2 Error 429: API rate limit exceeded during get {entrez_id} {variant} info, try again: {response.text}" + bcolors.ENDC)
                    time.sleep(random.uniform(0.25, 0.5))
                else:
                    print(bcolors.FAIL + f"Error {response.status_code}: {response.text}" + bcolors.ENDC)
                    time.sleep(random.uniform(0.25, 0.5))

    @staticmethod
    # Get gene information
    def all_variant(entrez_id, from_id=False):
        global headers

        if not from_id:
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

                all_variants = []

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

                if not from_id:
                    for elem in root.iter('Gene-commentary_accession'):
                        if elem.text.startswith('NC_'):
                            chromosome = elem.text
                            break

                    for variant in variants:
                        exon_coords = []  # Liste pour stocker toutes les paires (start, end)
                        found_variant = False
                        k_found = False
                        orientation = ""

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
                                start = int(elem.text)  # Ajuster pour 1-based
                            elif k_found and found_variant and elem.tag == "Seq-interval_to":
                                end = int(elem.text)  # Ajuster pour 1-based
                                exon_coords.append((start, end))  # Ajouter le couple (start, end) pour chaque exon
                            elif k_found and found_variant and elem.tag == "Na-strand" and orientation == "":
                                orientation += elem.attrib.get("value")

                            elif elem.tag == "Org-ref_taxname":
                                species_API = elem.text

                            elif elem.tag == 'Gene-ref_locus':
                                gene_name = elem.text

                        # Si des exons sont trouvés, normalisez les coordonnées
                        if exon_coords:

                            if orientation == "minus":
                                # Inverser les exons pour le brin négatif
                                exon_coords = [(end, start) for start, end in exon_coords]

                            # La première coordonnée du premier exon (start) devient la référence pour la normalisation
                            first_exon_start = exon_coords[0][0]

                            # Normalisation: soustraction de la première coordonnée
                            normalized_exon_coords = [
                                (abs(start - first_exon_start), abs(end - first_exon_start)) for start, end in
                                exon_coords]

                            # Ajouter les coordonnées chromosomiques et normalisées à la liste
                            all_variants.append(
                                (variant, gene_name, chraccver, exon_coords, normalized_exon_coords, species_API))

                    if len(all_variants) > 0:
                        print(
                            bcolors.OKGREEN + f"Response 200: Transcript(s) found(s) for {entrez_id}: {all_variants}" + bcolors.ENDC)
                        return all_variants, f"Transcript(s) found(s) for {entrez_id}: {all_variants}"
                    else:
                        all_variants.append(("Error 200", None, None, None, None, None))
                        print(
                            bcolors.WARNING + f"Error 200: Transcript not found(s) for {entrez_id}." + bcolors.ENDC)
                        return "Error 200", f"Transcript not found(s) for {entrez_id}."

                elif from_id:
                    if len(tv) > 0:
                        if "transcript variant 1" in tv:
                            associations = dict(zip(tv, variants))
                            return associations["transcript variant 1"]
                        else:
                            return variants[0]
                    elif len(tv) == 0 and len(variants) > 0:
                        return variants[0]
                    else:
                        return None

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