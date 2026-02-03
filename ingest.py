import json
import sys
from neo4j import GraphDatabase
from tqdm import tqdm
from rdkit import Chem, RDLogger
import os

# Suppress RDKit Warnings about hydrogen removal
# "WARNING: not removing hydrogen atom without neighbors"
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

# ----------------------------
URI = "bolt://localhost:7687"
USER = "neo4j"
PASSWORD = "chem_react"
DATA_DIR = "data"  # directory containing JSONL files
# ----------------------------

def canonicalize_smiles(smiles):
    """Convert SMILES to canonical form using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None

def get_smiles_from_identifiers(identifiers):
    for ident in identifiers:
        if ident.get("type") == "SMILES":
            return ident.get("value")
    return None

def extract_reactants(inputs):
    reactants = set()
    for group in inputs.values():
        for comp in group.get("components", []):
            if comp.get("reaction_role") != "REACTANT":
                continue
            smi = get_smiles_from_identifiers(comp.get("identifiers", []))
            if smi:
                canon = canonicalize_smiles(smi)
                if canon:
                    reactants.add(canon)
    return list(reactants)

def extract_products(outcomes):
    products = set()
    for outcome in outcomes:
        for prod in outcome.get("products", []):
            smi = get_smiles_from_identifiers(prod.get("identifiers", []))
            if smi:
                canon = canonicalize_smiles(smi)
                if canon:
                    products.add(canon)
    return list(products)

def extract_reaction_type(rxn):
    for ident in rxn.get("identifiers", []):
        if ident.get("type") == "REACTION_TYPE":
            return ident.get("value")
    return None


def extract_temperature(conditions):
    temp = conditions.get("temperature", {})
    setpoint = temp.get("setpoint")
    if setpoint:
        return (
            setpoint.get("value"),
            setpoint.get("units")
        )
    return (None, None)


def extract_solvents(inputs):
    solvents = set()
    for group in inputs.values():
        for comp in group.get("components", []):
            if comp.get("reaction_role") != "SOLVENT":
                continue
            smi = get_smiles_from_identifiers(comp.get("identifiers", []))
            if smi:
                canon = canonicalize_smiles(smi)
                if canon:
                    solvents.add(canon)
    return list(solvents)


def extract_catalysts(inputs):
    cats = set()
    for group in inputs.values():
        for comp in group.get("components", []):
            if comp.get("reaction_role") != "CATALYST":
                continue
            smi = get_smiles_from_identifiers(comp.get("identifiers", []))
            if smi:
                canon = canonicalize_smiles(smi)
                if canon:
                    cats.add(canon)
    return list(cats)


def extract_yield(outcomes):
    for outcome in outcomes:
        for prod in outcome.get("products", []):
            for meas in prod.get("measurements", []):
                if meas.get("type") == "YIELD":
                    pct = meas.get("percentage", {})
                    return pct.get("value")
    return None


def extract_procedure(rxn):
    return rxn.get("notes", {}).get("procedure_details")


def extract_publication_url(rxn):
    return rxn.get("provenance", {}).get("publication_url")

# ----------------------------

def ingest_batch(tx, batch):
    """
    Ingests a batch of reactions efficiently using UNWIND and CALL subqueries.
    """
    # Using CALL subqueries to isolate the cardinality of reactants vs products 
    # and prevent Cartesian explosions within the batch.
    query = """
    UNWIND $batch AS data
    MERGE (r:Reaction {reaction_id: data.rid})
    SET r.reaction_type = data.reaction_type,
        r.solvents = data.solvents,
        r.catalysts = data.catalysts,
        r.yield_percent = data.yield_percent,
        r.procedure = data.procedure
    
    WITH r, data
    CALL {
        WITH r, data
        UNWIND data.reactants AS rsmi
        MERGE (m1:Molecule {smiles: rsmi})
        MERGE (m1)-[:REACTANT_IN]->(r)
    }
    CALL {
        WITH r, data
        UNWIND data.products AS psmi
        MERGE (m2:Molecule {smiles: psmi})
        MERGE (r)-[:PRODUCES]->(m2)
    }
    """
    tx.run(query, batch=batch)


# ----------------------------

def main():
    # Find all JSONL files in DATA_DIR
    if not os.path.isdir(DATA_DIR):
        print(f"Data directory not found: {DATA_DIR}")
        sys.exit(1)

    jsonl_files = sorted([f for f in os.listdir(DATA_DIR) if f.endswith(".jsonl")])

    if not jsonl_files:
        print(f"No JSONL files found in {DATA_DIR}")
        sys.exit(1)

    print(f"Found {len(jsonl_files)} JSONL file(s) to ingest")

    driver = GraphDatabase.driver(URI, auth=(USER, PASSWORD))
    batch_size = 500

    try:
        for jsonl_filename in jsonl_files:
            jsonl_file = os.path.join(DATA_DIR, jsonl_filename)
            print(f"\nIngesting {jsonl_filename}...")

            batch = []
            with driver.session() as session:
                with open(jsonl_file) as f:
                    for line in tqdm(f):
                        rxn = json.loads(line)

                        # Get reaction ID, fallback if missing
                        rid = rxn.get("reaction_id")
                        if not rid:
                            for ident in rxn.get("identifiers", []):
                                if ident.get("type") == "CUSTOM":
                                    rid = ident.get("value")
                                    break
                        if not rid:
                            continue

                        inputs = rxn.get("inputs", {})
                        outcomes = rxn.get("outcomes", [])

                        batch.append({
                            "rid": rid,
                            "reactants": extract_reactants(inputs),
                            "products": extract_products(outcomes),
                            "reaction_type": extract_reaction_type(rxn),
                            "solvents": extract_solvents(inputs),
                            "catalysts": extract_catalysts(inputs),
                            "yield_percent": extract_yield(outcomes),
                            "procedure": extract_procedure(rxn)
                        })

                        if len(batch) >= batch_size:
                            session.execute_write(ingest_batch, batch)
                            batch = []

                    # Final partial batch
                    if batch:
                        session.execute_write(ingest_batch, batch)

            print(f"Completed {jsonl_filename}")
    finally:
        driver.close()
        print("\nIngestion complete.")

# ----------------------------

if __name__ == "__main__":
    main()
