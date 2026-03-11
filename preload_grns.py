# preload_grns.py

import os
import celloracle as co

# Correct environment variable
CACHE_DIR = os.path.join(os.getcwd(), "celloracle_cache")
os.environ["CELLORACLE_DATA_DIR"] = CACHE_DIR
os.makedirs(os.path.join(CACHE_DIR, "promoter_base_GRN"), exist_ok=True)

# Species to GRN loader mapping
species_to_loader = {
    "human": co.data.load_human_promoter_base_GRN,
    "mouse": co.data.load_mouse_promoter_base_GRN,
    "pig": co.data.load_Pig_promoter_base_GRN,
    "rat": co.data.load_rat_promoter_base_GRN,
    "zebrafish": co.data.load_zebrafish_promoter_base_GRN,
    "chicken": co.data.load_chicken_promoter_base_GRN,
    "celegans": co.data.load_Celegans_promoter_base_GRN,
    "drosophila": co.data.load_drosophila_promoter_base_GRN,
    "scerevisiae": co.data.load_Scerevisiae_promoter_base_GRN,
    "arabidopsis": co.data.load_arabidopsis_promoter_base_GRN,
    "xenopus_laevis": co.data.load_xenopus_laevis_promoter_base_GRN,
    "xenopus_tropicalis": co.data.load_xenopus_tropicalis_promoter_base_GRN,
}

for species, loader_fn in species_to_loader.items():
    print(f"🔽 Downloading base GRN for {species}...")
    try:
        loader_fn()
        print(f"✅ {species}: Downloaded and cached.")
    except Exception as e:
        print(f"❌ Failed to load {species}: {e}")
