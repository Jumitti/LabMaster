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


ucsc_species = {
    "Homo sapiens": {'org': 'Human', 'db': 'hg38', 'wp_target': ['genome', 'hg38KgSeqV47']},
    "Mus musculus": {'org': 'Mouse', 'db': 'mm39', 'wp_target': ['genome', 'mm39KgSeqVM36']},
    "Anopheles gambiae": {'org': 'A.+gambiae', 'db': 'anoGam3', 'wp_target': ['genome']},
    "Apis mellifera": {'org': 'A.+mellifera', 'db': 'apiMel2', 'wp_target': ['genome']},
    "Xenopus laevis": {'org': 'African+clawed+frog', 'db': 'xenLae2', 'wp_target': ['genome']},
    "Vicugna pacos": {'org': 'Alpaca', 'db': 'vicPac2', 'wp_target': ['genome']},
    "Alligator mississippiensis": {'org': 'American+alligator', 'db': 'allMis1', 'wp_target': ['genome']},
    "Dasypus novemcinctus": {'org': 'Armadillo', 'db': 'dasNov3', 'wp_target': ['genome']},
    "Gadus morhua": {'org': 'Atlantic+cod', 'db': 'gadMor1', 'wp_target': ['genome']},
    "Papio anubis": {'org': 'Baboon', 'db': 'papAnu4', 'wp_target': ['genome']},
    "Bison bison": {'org': 'Bison', 'db': 'bisBis1', 'wp_target': ['genome']},
    "Pan paniscus": {'org': 'Bonobo', 'db': 'panPan3', 'wp_target': ['genome']},
    "Apteryx australis": {'org': 'Brown+kiwi', 'db': 'aptMan1', 'wp_target': ['genome']},
    "Melopsittacus undulatus": {'org': 'Budgerigar', 'db': 'melUnd1', 'wp_target': ['genome']},
    "Otolemur garnettii": {'org': 'Bushbaby', 'db': 'otoGar3', 'wp_target': ['genome']},
    "Caenorhabditis brenneri": {'org': 'C.+brenneri', 'db': 'caePb2', 'wp_target': ['genome']},
    "Caenorhabditis briggsae": {'org': 'C.+briggsae', 'db': 'cb3', 'wp_target': ['genome']},
    "Caenorhabditis elegans": {'org': 'C.+elegans', 'db': 'ce11', 'wp_target': ['genome']},
    "Ciona intestinalis": {'org': 'C.+intestinalis', 'db': 'ci3', 'wp_target': ['genome']},
    "Caenorhabditis japonica": {'org': 'C.+japonica', 'db': 'caeJap1', 'wp_target': ['genome']},
    "Caenorhabditis remanei": {'org': 'C.+remanei', 'db': 'caeRem3', 'wp_target': ['genome']},
    "Felis catus": {'org': 'Cat', 'db': 'felCat9', 'wp_target': ['genome']},
    "Gallus gallus": {'org': 'Chicken', 'db': 'galGal6', 'wp_target': ['genome']},
    "Pan troglodytes": {'org': 'Chimp', 'db': 'panTro6', 'wp_target': ['genome']},
    "Cricetulus griseus": {'org': 'Chinese+hamster', 'db': 'criGriChoV2', 'wp_target': ['genome']},
    "Manis pentadactyla": {'org': 'Chinese+pangolin', 'db': 'manPen1', 'wp_target': ['genome']},
    "Latimeria chalumnae": {'org': 'Coelacanth', 'db': 'latCha1', 'wp_target': ['genome']},
    "Bos taurus": {'org': 'Cow', 'db': 'bosTau9', 'wp_target': ['genome']},
    "Macaca fascicularis": {'org': 'Crab-eating+macaque', 'db': 'macFas5', 'wp_target': ['genome']},
    "Drosophila ananassae": {'org': 'D.+ananassae', 'db': 'droAna2', 'wp_target': ['genome']},
    "Drosophila erecta": {'org': 'D.+erecta', 'db': 'droEre1', 'wp_target': ['genome']},
    "Drosophila grimshawi": {'org': 'D.+grimshawi', 'db': 'droGri1', 'wp_target': ['genome']},
    "Drosophila melanogaster": {'org': 'D.+melanogaster', 'db': 'dm6', 'wp_target': ['genome']},
    "Drosophila mojavensis": {'org': 'D.+mojavensis', 'db': 'droMoj2', 'wp_target': ['genome']},
    "Drosophila persimilis": {'org': 'D.+persimilis', 'db': 'droPer1', 'wp_target': ['genome']},
    "Drosophila pseudoobscura": {'org': 'D.+pseudoobscura', 'db': 'dp3', 'wp_target': ['genome']},
    "Drosophila sechellia": {'org': 'D.+sechellia', 'db': 'droSec1', 'wp_target': ['genome']},
    "Drosophila simulans": {'org': 'D.+simulans', 'db': 'droSim1', 'wp_target': ['genome']},
    "Drosophila virilis": {'org': 'D.+virilis', 'db': 'droVir2', 'wp_target': ['genome']},
    "Drosophila yakuba": {'org': 'D.+yakuba', 'db': 'droYak2', 'wp_target': ['genome']},
    "Canis lupus familiaris": {'org': 'Dog', 'db': 'canFam6', 'wp_target': ['genome']},
    "Delphinidae": {'org': 'Dolphin', 'db': 'turTru2', 'wp_target': ['genome']},
    "Zaire ebolavirus": {'org': 'Ebola+virus', 'db': 'eboVir3', 'wp_target': ['genome']},
    "Loxodonta africana": {'org': 'Elephant', 'db': 'loxAfr3', 'wp_target': ['genome']},
    "Callorhinchus milii": {'org': 'Elephant+shark', 'db': 'calMil1', 'wp_target': ['genome']},
    "Mustela putorius furo": {'org': 'Ferret+', 'db': 'musFur1', 'wp_target': ['genome']},
    "Takifugu rubripes": {'org': 'Fugu', 'db': 'fr3', 'wp_target': ['genome']},
    "Thamnophis sirtalis": {'org': 'Garter+snake', 'db': 'thaSir1', 'wp_target': ['genome']},
    "Hylobates": {'org': 'Gibbon', 'db': 'nomLeu3', 'wp_target': ['genome']},
    "Aquila chrysaetos": {'org': 'Golden+eagle', 'db': 'aquChr2', 'wp_target': ['genome']},
    "Rhinopithecus roxellana": {'org': 'Golden+snub-nosed+monkey', 'db': 'rhiRox1', 'wp_target': ['genome']},
    "Gorilla gorilla": {'org': 'Gorilla', 'db': 'gorGor6', 'wp_target': ['genome']},
    "Chlorocebus sabaeus": {'org': 'Green+monkey', 'db': 'chlSab2', 'wp_target': ['genome']},
    "Cavia porcellus": {'org': 'Guinea+pig', 'db': 'cavPor3', 'wp_target': ['genome']},
    "Neomonachus schauinslandi": {'org': 'Hawaiian+monk+seal', 'db': 'neoSch1', 'wp_target': ['genome']},
    "Erinaceus europaeus": {'org': 'Hedgehog', 'db': 'eriEur2', 'wp_target': ['genome']},
    "Equus caballus": {'org': 'Horse', 'db': 'equCab3', 'wp_target': ['genome']},
    "Dipodomys": {'org': 'Kangaroo+rat', 'db': 'dipOrd1', 'wp_target': ['genome']},
    "Petromyzon marinus": {'org': 'Lamprey', 'db': 'petMar2', 'wp_target': ['genome']},
    "Branchiostoma": {'org': 'Lancelet', 'db': 'braFlo1', 'wp_target': ['genome']},
    "Sceloporus": {'org': 'Lizard', 'db': 'anoCar2', 'wp_target': ['genome']},
    "Cynocephalus variegatus": {'org': 'Malayan+flying+lemur', 'db': 'galVar1', 'wp_target': ['genome']},
    "Trichechus manatus": {'org': 'Manatee', 'db': 'triMan1', 'wp_target': ['genome']},
    "Callithrix jacchus": {'org': 'Marmoset', 'db': 'calJac4', 'wp_target': ['genome']},
    "Oryzias latipes": {'org': 'Medaka', 'db': 'oryLat2', 'wp_target': ['genome']},
    "Geospiza fortis": {'org': 'Medium+ground+finch', 'db': 'geoFor1', 'wp_target': ['genome']},
    "Pteropus": {'org': 'Megabat', 'db': 'pteVam1', 'wp_target': ['genome']},
    "Myotis lucifugus": {'org': 'Little+brown+bat', 'db': 'myoLuc2', 'wp_target': ['genome']},
    "Balaenoptera acutorostrata": {'org': 'Minke+whale', 'db': 'balAcu1', 'wp_target': ['genome']},
    "Microcebus murinus": {'org': 'Mouse+lemur', 'db': 'micMur2', 'wp_target': ['genome']},
    "Heterocephalus glaber": {'org': 'Naked+mole-rat', 'db': 'hetGla2', 'wp_target': ['genome']},
    "Oreochromis niloticus": {'org': 'Nile+tilapia', 'db': 'oreNil2', 'wp_target': ['genome']},
    "Monodelphis domestica": {'org': 'Opossum', 'db': 'monDom5', 'wp_target': ['genome']},
    "Pongo pygmaeus": {'org': 'Orangutan', 'db': 'ponAbe3', 'wp_target': ['genome']},
    "Pristionchus pacificus": {'org': 'P.+pacificus', 'db': 'priPac1', 'wp_target': ['genome']},
    "Chrysemys picta": {'org': 'Painted+turtle', 'db': 'chrPic1', 'wp_target': ['genome']},
    "Ailuropoda melanoleuca": {'org': 'Panda', 'db': 'ailMel1', 'wp_target': ['genome']},
    "Sus scrofa": {'org': 'Pig', 'db': 'susScr11', 'wp_target': ['genome']},
    "Ochotona princeps": {'org': 'Pika', 'db': 'ochPri3', 'wp_target': ['genome']},
    "Ornithorhynchus anatinus": {'org': 'Platypus', 'db': 'ornAna2', 'wp_target': ['genome']},
    "Nasalis larvatus": {'org': 'Proboscis+monkey', 'db': 'nasLar1', 'wp_target': ['genome']},
    "Oryctolagus cuniculus": {'org': 'Rabbit', 'db': 'oryCun2', 'wp_target': ['genome']},
    "Rattus norvegicus": {'org': 'Rat', 'db': 'rn7', 'wp_target': ['genome']},
    "Macaca mulatta": {'org': 'Rhesus', 'db': 'rheMac10', 'wp_target': ['genome']},
    "Procavia capensis": {'org': 'Rock+hyrax', 'db': 'proCap1', 'wp_target': ['genome']},
    "Saccharomyces cerevisiae": {'org': 'S.+cerevisiae', 'db': 'sacCer3', 'wp_target': ['genome']},
    "Strongylocentrotus purpuratus": {'org': 'S.+purpuratus', 'db': 'strPur2', 'wp_target': ['genome']},
    "Aplysia californica": {'org': 'Sea+hare', 'db': 'aplCal1', 'wp_target': ['genome']},
    "Enhydra lutris nereis": {'org': 'Southern+sea+otter', 'db': 'enhLutNer1', 'wp_target': ['genome']},
    "Ovis aries": {'org': 'Sheep', 'db': 'oviAri4', 'wp_target': ['genome']},
    "Sorex araneus": {'org': 'Shrew', 'db': 'sorAra2', 'wp_target': ['genome']},
    "Bradypus": {'org': 'Sloth', 'db': 'choHof1', 'wp_target': ['genome']},
    "Sciurus": {'org': 'Squirrel', 'db': 'speTri2', 'wp_target': ['genome']},
    "Saimiri sciureus": {'org': 'Squirrel+monkey', 'db': 'saiBol1', 'wp_target': ['genome']},
    "Gasterosteus aculeatus": {'org': 'Stickleback', 'db': 'gasAcu1', 'wp_target': ['genome']},
    "Tarsius syrichta": {'org': 'Tarsier', 'db': 'tarSyr2', 'wp_target': ['genome']},
    "Sarcophilus harrisii": {'org': 'Tasmanian+devil', 'db': 'sarHar1', 'wp_target': ['genome']},
    "Tenrec ecaudatus": {'org': 'Tenrec', 'db': 'echTel2', 'wp_target': ['genome']},
    "Tetraodon nigroviridis": {'org': 'Tetraodon', 'db': 'tetNig2', 'wp_target': ['genome']},
    "Nanorana parkeri": {'org': 'Tibetan+frog', 'db': 'nanPar1', 'wp_target': ['genome']},
    "Tupaia": {'org': 'Tree+shrew', 'db': 'tupBel1', 'wp_target': ['genome']},
    "Meleagris gallopavo": {'org': 'Turkey', 'db': 'melGal5', 'wp_target': ['genome']},
    "Severe acute respiratory syndrome coronavirus 2": {'org': 'SARS-CoV-2', 'db': 'wuhCor1', 'wp_target': ['genome']},
    "Monkeypox virus": {'org': 'Monkeypox+virus', 'db': 'mpxvRivers', 'wp_target': ['genome']},
    "Macropus": {'org': 'Wallaby', 'db': 'macEug2', 'wp_target': ['genome']},
    "Ceratotherium simum": {'org': 'White+rhinoceros', 'db': 'cerSim1', 'wp_target': ['genome']},
    "Xenopus tropicalis": {'org': 'X.+tropicalis', 'db': 'xenTro10', 'wp_target': ['genome']},
    "Taeniopygia guttata": {'org': 'Zebra+finch', 'db': 'taeGut2', 'wp_target': ['genome']},
    "Danio rerio": {'org': 'Zebrafish', 'db': 'danRer11', 'wp_target': ['genome']},
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
    def design_primers(variant, gene_name, species, sequence, exons, nb_primers, product_size_settings, ucsc_validation=False, only_validated=False):
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
                'PRIMER_MIN_TM': 57.0,  # Minimum melting temperature (°C)
                'PRIMER_MAX_TM': 63.0,  # Maximum melting temperature (°C)
                'PRIMER_MIN_GC': 40.0,  # Minimum GC percentage (%)
                'PRIMER_MAX_GC': 60.0,  # Maximum GC percentage (%)
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
                'PRIMER_PRODUCT_SIZE_RANGE': [product_size_settings],  # Product size range

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

                                        if species in ucsc_species.keys() and ucsc_validation is True:
                                            org = ucsc_species[species]["org"]
                                            db = ucsc_species[species]["db"]
                                            wp_targets = ucsc_species[species]["wp_target"]
                                            validation_relative, validation_absolute = NCBIdna.fetch_ucsc_pcr_results(
                                                species, org, db, wp_targets, left_seq, right_seq,
                                                amplicon_size_abs, product_size_settings[1])
                                        else:
                                            validation_relative, validation_absolute = None, None

                                        if only_validated is True:
                                            if not validation_relative and not validation_absolute:
                                                continue

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
                                            'validation_relative': validation_relative,
                                            'validation_absolute': validation_absolute,
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

    @staticmethod
    def fetch_ucsc_pcr_results(species, org, db, wp_targets, wp_f, wp_r, amplicon_size_abs, max_product_size):
        validation_relative = None
        validation_absolute = None

        for wp_target in wp_targets:
            if wp_target == "genome":
                amplicon_size = amplicon_size_abs + 100
            else:
                amplicon_size = 100 + max_product_size

            base_url = "https://genome.ucsc.edu/cgi-bin/hgPcr"
            params = {
                "org": org,
                "db": db,
                "wp_target": wp_target,
                "wp_f": wp_f,
                "wp_r": wp_r,
                "Submit": "Submit",
                "wp_size": amplicon_size * 2,
                "wp_perfect": 15,
                "wp_good": 15,
                "boolshad.wp_flipReverse": 0,
                "wp_append": "on",
                "boolshad.wp_append": 0,
            }

            response = requests.get(base_url, params=params)

            if response.status_code != 200:
                if wp_target == "genome":
                    validation_absolute = f"Error {response.status_code}"
                else:
                    validation_relative = f"Error {response.status_code}"

            soup = BeautifulSoup(response.text, "html.parser")
            result_section = soup.find("pre")

            if result_section:
                results = result_section.text.strip().split("\n")
                fragment_sizes = []

                for line in results:
                    match = re.search(r'(\d+)bp', line)
                    if match:
                        size = int(match.group(1))
                        if size <= amplicon_size:
                            fragment_sizes.append(size)

                if not fragment_sizes:
                    if wp_target == "genome":
                        validation_absolute = "Not found"
                    else:
                        validation_relative = "Not found"

                if all(size == fragment_sizes[0] for size in fragment_sizes):
                    if wp_target == "genome":
                        validation_absolute = True
                    else:
                        validation_relative = True
                else:
                    if wp_target == "genome":
                        validation_absolute = False
                    else:
                        validation_relative = False
            else:
                if wp_target == "genome":
                    validation_absolute = "Not found"
                else:
                    validation_relative = "Not found"

        if validation_relative is None:
            ncbi_pcr_results = NCBIdna.ncb_pcr_in_silico(species, wp_f, wp_r)

            filtered_lengths = [res["product_length"] for res in ncbi_pcr_results if
                                res["product_length"] is not None and res["product_length"] < 100 + max_product_size]

            if len(set(filtered_lengths)) <= 1:
                validation_relative = True
            else:
                validation_relative = False

        return validation_relative, validation_absolute

    @staticmethod
    def ncb_pcr_in_silico(species, primer_fwd, primer_rev):
        url = (
            "https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?"
            "CMD=request&CON_ANEAL_OLIGO=50.0&CON_DNTPS=0.6&DIVA_CATIONS=1.5&EVALUE=30000"
            "&GC_CLAMP=0&HITSIZE=50000&LOW_COMPLEXITY_FILTER=on&MAX_CANDIDATE_PRIMER=500"
            "&MAX_INTRON_SIZE=1000000&MAX_TARGET_PER_TEMPLATE=100&MAX_TARGET_SIZE=4000"
            "&MIN_INTRON_SIZE=1000&MISMATCH_REGION_LENGTH=5&MONO_CATIONS=50.0"
            f"&ORGANISM={species}&OVERLAP_3END=4&OVERLAP_5END=7&POLYX=5"
            "&PRIMER_3END_SPECIFICITY_MISMATCH=1&PRIMER_INTERNAL_OLIGO_MAX_GC=80.0"
            "&PRIMER_INTERNAL_OLIGO_MAX_SIZE=27&PRIMER_INTERNAL_OLIGO_MAX_TM=63.0"
            "&PRIMER_INTERNAL_OLIGO_MIN_GC=20.0&PRIMER_INTERNAL_OLIGO_MIN_SIZE=18"
            "&PRIMER_INTERNAL_OLIGO_MIN_TM=57.0&PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT=50"
            "&PRIMER_INTERNAL_OLIGO_OPT_SIZE=20&PRIMER_INTERNAL_OLIGO_OPT_TM=60.0"
            f"&PRIMER_LEFT_INPUT={primer_fwd}&PRIMER_MAX_DIFF_TM=3&PRIMER_MAX_END_GC=5"
            "&PRIMER_MAX_END_STABILITY=9&PRIMER_MAX_GC=80.0&PRIMER_MAX_HAIRPIN_TH=24.0"
            "&PRIMER_MAX_SELF_ANY_TH=45.0&PRIMER_MAX_SELF_END_TH=35.0"
            "&PRIMER_MAX_SIZE=25&PRIMER_MAX_TEMPLATE_MISPRIMING=12.00"
            "&PRIMER_MAX_TEMPLATE_MISPRIMING_TH=40.00&PRIMER_MAX_TM=63.0"
            "&PRIMER_MIN_GC=20.0&PRIMER_MIN_SIZE=15&PRIMER_MIN_TM=57.0"
            "&PRIMER_MISPRIMING_LIBRARY=AUTO&PRIMER_NUM_RETURN=10&PRIMER_ON_SPLICE_SITE=0"
            "&PRIMER_OPT_SIZE=20&PRIMER_OPT_TM=60.0&PRIMER_PAIR_MAX_COMPL_ANY=8.00"
            "&PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0&PRIMER_PAIR_MAX_COMPL_END=3.00"
            "&PRIMER_PAIR_MAX_COMPL_END_TH=35.0&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00"
            "&PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=70.00&PRIMER_PRODUCT_MAX=1000"
            "&PRIMER_PRODUCT_MIN=70"
            f"&PRIMER_RIGHT_INPUT={primer_rev}&PRIMER_SPECIFICITY_DATABASE=refseq_mrna"
            "&SALT_FORMULAR=1&SEARCH_SPECIFIC_PRIMER=on&SEARCHMODE=0&SELF_ANY=8.00"
            "&SELF_END=3.00&SHOW_SVIEWER=on&SPLICE_SITE_OVERLAP_3END=4"
            "&SPLICE_SITE_OVERLAP_3END_MAX=8&SPLICE_SITE_OVERLAP_5END=7&TM_METHOD=1"
            "&TOTAL_MISMATCH_IGNORE=6&TOTAL_PRIMER_SPECIFICITY_MISMATCH=1&UNGAPPED_BLAST=on"
            "&USER_TYPE=2&WORD_SIZE=7"
        )

        session = requests.Session()
        response = session.get(url)

        if response.status_code != 200:
            print(f"❌ Erreur lors de la requête : {response.status_code}")
            exit()

        soup = BeautifulSoup(response.text, "html.parser")
        job_id_input = soup.find("input", {"name": "job_key"})
        if job_id_input:
            job_key = job_id_input.get("value")
        else:
            match = re.search(r"job_key=([A-Za-z0-9_-]+)", response.text)
            if match:
                job_key = match.group(1)
            else:
                print("❌ Impossible de récupérer le JOB ID.")
                exit()

        print(f"✅ JOB ID récupéré : {job_key}")

        status_url = f"https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key={job_key}&CMD=get"

        while True:
            time.sleep(10)
            status_response = session.get(status_url)

            soup = BeautifulSoup(status_response.text, "html.parser")
            print(soup)

            primer_pair = soup.find("a", {"name": "0"})
            forward_primer = soup.find("th", string="Forward primer")
            reverse_primer = soup.find("th", string="Reverse primer")

            if primer_pair and forward_primer and reverse_primer:
                print("✅ Résultats prêts !")
                break
            else:
                print("⏳ En attente des résultats...")

        results = []

        entries = soup.find_all("a", href=re.compile(r"viewer.fcgi\?db=nucleotide"))

        for entry in entries:
            gene_name = entry.text.strip()
            gene_info = entry.next_sibling.strip() if entry.next_sibling else "N/A"

            pre_block = entry.find_next("pre").text.strip()

            product_length_match = re.search(r"product length = (\d+)", pre_block)
            product_length = int(product_length_match.group(1)) if product_length_match else None

            forward_primer, forward_start, forward_template, forward_end = "N/A", "N/A", "N/A", "N/A"
            reverse_primer, reverse_start, reverse_template, reverse_end = "N/A", "N/A", "N/A", "N/A"

            forward_match = re.search(r"Forward primer\s+\d+\s+([A-Za-z]+)\s+\d+", pre_block)
            reverse_match = re.search(r"Reverse primer\s+\d+\s+([A-Za-z]+)\s+\d+", pre_block)

            if forward_match:
                forward_primer = forward_match.group(1)
            if reverse_match:
                reverse_primer = reverse_match.group(1)

            template_matches = re.findall(r"Template\s+(\d+)\s+([A-Za-z.\s]+)\s+(\d+)", pre_block)

            if len(template_matches) > 0:
                forward_start, forward_template, forward_end = template_matches[0]
            if len(template_matches) > 1:
                reverse_start, reverse_template, reverse_end = template_matches[1]

            results.append({
                "gene": gene_info,
                "name": gene_name,
                "product_length": product_length,
                "forward_primer": forward_primer,
                "forward_start": forward_start,
                "forward_template": forward_template,
                "forward_end": forward_end,
                "reverse_primer": reverse_primer,
                "reverse_start": reverse_start,
                "reverse_template": reverse_template,
                "reverse_end": reverse_end,
            })

        return results
