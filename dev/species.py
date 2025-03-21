import re

# Liste des espèces à traiter
species_list = [
    "Homo sapiens",
    "Mus musculus",
    "Anopheles gambiae",
    "Apis mellifera",
    "Xenopus laevis",
    "Vicugna pacos",
    "Alligator mississippiensis",
    "Dasypus novemcinctus",
    "Gadus morhua",
    "Papio anubis",
    "Bison bison",
    "Pan paniscus",
    "Apteryx australis",
    "Melopsittacus undulatus",
    "Otolemur garnettii",
    "Caenorhabditis brenneri",
    "Caenorhabditis briggsae",
    "Caenorhabditis elegans",
    "Ciona intestinalis",
    "Caenorhabditis japonica",
    "Caenorhabditis remanei",
    "Felis catus",
    "Gallus gallus",
    "Pan troglodytes",
    "Cricetulus griseus",
    "Manis pentadactyla",
    "Latimeria chalumnae",
    "Bos taurus",
    "Macaca fascicularis",
    "Drosophila ananassae",
    "Drosophila erecta",
    "Drosophila grimshawi",
    "Drosophila melanogaster",
    "Drosophila mojavensis",
    "Drosophila persimilis",
    "Drosophila pseudoobscura",
    "Drosophila sechellia",
    "Drosophila simulans",
    "Drosophila virilis",
    "Drosophila yakuba",
    "Canis lupus familiaris",
    "Delphinidae",
    "Zaire ebolavirus",
    "Loxodonta africana",
    "Callorhinchus milii",
    "Mustela putorius furo",
    "Takifugu rubripes",
    "Thamnophis sirtalis",
    "Hylobates",
    "Aquila chrysaetos",
    "Rhinopithecus roxellana",
    "Gorilla gorilla",
    "Chlorocebus sabaeus",
    "Cavia porcellus",
    "Neomonachus schauinslandi",
    "Erinaceus europaeus",
    "Equus caballus",
    "Dipodomys",
    "Petromyzon marinus",
    "Branchiostoma",
    "Sceloporus",
    "Cynocephalus variegatus",
    "Trichechus manatus",
    "Callithrix jacchus",
    "Oryzias latipes",
    "Geospiza fortis",
    "Pteropus",
    "Myotis lucifugus",
    "Balaenoptera acutorostrata",
    "Microcebus murinus",
    "Heterocephalus glaber",
    "Oreochromis niloticus",
    "Monodelphis domestica",
    "Pongo pygmaeus",
    "Pristionchus pacificus",
    "Chrysemys picta",
    "Ailuropoda melanoleuca",
    "Sus scrofa",
    "Ochotona princeps",
    "Ornithorhynchus anatinus",
    "Nasalis larvatus",
    "Oryctolagus cuniculus",
    "Rattus norvegicus",
    "Macaca mulatta",
    "Procavia capensis",
    "Saccharomyces cerevisiae",
    "Strongylocentrotus purpuratus",
    "Aplysia californica",
    "Enhydra lutris nereis",
    "Ovis aries",
    "Sorex araneus",
    "Bradypus",
    "Sciurus",
    "Saimiri sciureus",
    "Gasterosteus aculeatus",
    "Tarsius syrichta",
    "Sarcophilus harrisii",
    "Tenrec ecaudatus",
    "Tetraodon nigroviridis",
    "Nanorana parkeri",
    "Tupaia",
    "Meleagris gallopavo",
    "Severe acute respiratory syndrome coronavirus 2",
    "Monkeypox virus",
    "Macropus",
    "Ceratotherium simum",
    "Xenopus tropicalis",
    "Taeniopygia guttata",
    "Danio rerio"
]

# Dictionnaire final
ucsc_species = {}

# Regex pour extraire les paramètres de l'URL
param_regex = re.compile(r'([\w\d_-]+)=([^&]+)')

for species in species_list:
    url = input(f"Entrez le lien UCSC pour {species}: ").strip()

    # Extraction des paramètres de l'URL
    params = dict(param_regex.findall(url))

    # Récupération des valeurs nécessaires
    org = params.get("org", "Unknown")
    db = params.get("db", "Unknown")
    wp_target = params.get("wp_target", "genome")

    # Vérification et formatage de wp_target
    if wp_target != "genome":
        wp_target = ["genome", wp_target]
    else:
        wp_target = ["genome"]

    # Stockage dans le dictionnaire final
    ucsc_species[species] = {"org": org, "db": db, "wp_target": wp_target}

# Affichage du dictionnaire final
print("\nucsc_species = {")
for species, values in ucsc_species.items():
    print(f'    "{species}": {values},')
print("}")
