import requests
from typing import List, Optional
from pydantic import BaseModel, Field
from rdkit import Chem, RDLogger
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate, FewShotChatMessagePromptTemplate
from langchain_core.runnables import RunnableLambda

load_dotenv()

#----------------------------
# Helpers for SMILES and Name Resolution
#----------------------------

def canonicalize_smiles(smiles: str) -> Optional[str]:
    """Validate and canonicalize a SMILES string, suppressing RDKit logs."""
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True) if mol else None
    except:
        return None
    finally:
        lg.setLevel(RDLogger.WARNING)

def resolve_name_to_smiles(name: str) -> Optional[str]:
    """Query OPSIN and PubChem to convert chemical name to SMILES."""
    # 1. Try OPSIN (Best for systematic names)
    try:
        r = requests.get(f"https://opsin.ch.cam.ac.uk/opsin/{name}.json", timeout=5)
        if r.status_code == 200:
            res = r.json().get("smiles")
            if res: return res
    except Exception:
        pass

    # 2. Try PubChem (Best for common names / drugs)
    try:
        # Request both Canonical and Isomeric SMILES
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,IsomericSMILES/JSON"
        r = requests.get(url, timeout=5)
        if r.status_code == 200:
            data = r.json()
            props = data.get("PropertyTable", {}).get("Properties", [])
            if props:
                # PubChem can return different keys. Prioritize Isomeric > Canonical > Connectivity > SMILES
                rec = props[0]
                return (
                    rec.get("IsomericSMILES") or 
                    rec.get("CanonicalSMILES") or 
                    rec.get("ConnectivitySMILES") or
                    rec.get("SMILES")
                )
    except Exception:
        pass
        
    return None

def resolve_identifier(identifier: str) -> Optional[str]:
    """Try to resolve input as a valid SMILES first, then as a Chemical Name."""
    smi = canonicalize_smiles(identifier)
    if smi: return smi
    
    smi_from_name = resolve_name_to_smiles(identifier)
    if smi_from_name:
        return canonicalize_smiles(smi_from_name)
    return None


#----------------------------
# Models & Agent
#----------------------------

class MoleculeEntity(BaseModel):
    smiles: Optional[str] = Field(None, description="The SMILES string if available")
    name: Optional[str] = Field(None, description="The common name if available")

    def resolve(self) -> Optional[str]:
        """Attempt to resolve this entity to a SMILES string."""
        # Prioritize name resolution as it leans on external APIs (Ground Truth)
        # The LLM often hallucinates SMILES for common names.
        if self.name:
            res = resolve_identifier(self.name)
            if res: return res
            
        if self.smiles:
            res = resolve_identifier(self.smiles)
            if res: return res
            
        return None

class SynthesisQuery(BaseModel):
    product: MoleculeEntity
    starting_materials: List[MoleculeEntity]

class ChemicalParserAgent:
    def __init__(self, model_name: str = "gpt-4o-mini", temperature: float = 0):
        self.llm = ChatOpenAI(model=model_name, temperature=temperature)
        self._setup_chain()

    def _setup_chain(self):
        # Few-shot prompting examples
        examples = [
            {"input": "How can I synthesize aspirin from salicylic acid and acetic anhydride?", 
             "output": {"product": {"name": "aspirin"}, "starting_materials": [{"name": "salicylic acid"}, {"name": "acetic anhydride"}]}},
            {"input": "Find a path to make C1=CC=CC=C1 (benzene) using acetylene", 
             "output": {"product": {"name": "benzene", "smiles": "C1=CC=CC=C1"}, "starting_materials": [{"name": "acetylene"}]}},
            {"input": "I have CCO and CC(=O)O. Can I make ethyl acetate from these?", 
             "output": {"product": {"name": "ethyl acetate"}, "starting_materials": [{"smiles": "CCO"}, {"smiles": "CC(=O)O"}]}}
        ]

        few_shot_prompt = FewShotChatMessagePromptTemplate(
            example_prompt=ChatPromptTemplate.from_messages([("human", "{input}"), ("ai", "{output}")]),
            examples=examples,
        )

        final_prompt = ChatPromptTemplate.from_messages([
            ("system", "You are an expert chemist. Extract the target product and starting reactants from the user's query."),
            few_shot_prompt,
            ("human", "{query}"),
        ])
        
        self.chain = final_prompt | self.llm.with_structured_output(SynthesisQuery)

    def get_runnable(self):
        """
        Returns a generic Runnable that accepts a query string (or dict with 'query')
        and returns { "target": str, "starting_materials": Set[str] }.
        """
        return RunnableLambda(self._runnable_entry)

    def _runnable_entry(self, input_data):
        # Handle both string input and dict input
        if isinstance(input_data, dict):
            query = input_data.get("query")
        else:
            query = str(input_data)
            
        return self.parse(query)

    def parse(self, query: str) -> dict:
        """
        Parses a natural language query using an LLM to extract target and starting materials.
        Returns a dict: { "target": ..., "starting_materials": ... } or None on failure.
        """
        try:
            result = self.chain.invoke({"query": query})
        except Exception as e:
            print(f"LLM Interaction Error: {e}")
            return None
        
        if not result:
            print("Error: Failed to parse query structure.")
            return None

        # Resolve product
        target_smiles = result.product.resolve()
        if not target_smiles:
            print(f"Error: Could not resolve target: {result.product}")
            return None

        # Resolve starting materials
        starting_materials = set()
        for m in result.starting_materials:
            smi = m.resolve()
            if smi:
                starting_materials.add(smi)
            else:
                print(f"Warning: Could not resolve starting material: {m}")
        
        if not starting_materials:
            print("Error: No valid starting molecules identified.")
            return None

        return {
            "target": target_smiles,
            "starting_materials": list(starting_materials) # Convert to list for easier serialization
        }

    target_molecule = "COc1cc2c(cc1Nc1ncc(C)c(-c3cnn4ccccc34)n1)N(C(C)=O)CC2"
    starting_molecules = set([
        "Cc1cnc(Cl)nc1-c1cnn2ccccc12",
        "CCO",
        "COc1cc2c(cc1[N+](=O)[O-])N(C(C)=O)CC2"
    ])
