from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnableLambda
from langchain_core.output_parsers import StrOutputParser
from dotenv import load_dotenv
import json

load_dotenv()

class ChemicalNarratorAgent:
    def __init__(self, model_name: str = "gpt-4o-mini", temperature: float = 0):
        self.llm = ChatOpenAI(model=model_name, temperature=temperature)
        self._setup_chain()

    def _setup_chain(self):
        template = """You are an expert chemist explaining a synthesis procedure.
You are given a list of sequential reaction steps for synthesizing a target molecule.
Some steps might have detailed procedure text from a database, others might be sparse.

Your goal is to write a cohesive, step-by-step "Lab Manual" style guide for this synthesis.
- Use imperative mood (e.g., "Add X to Y", "Heat to 50C").
- If procedure text is available, paraphrase it clearly.
- If procedure text is missing ("No procedure details found."), infer standard techniques for the reaction type if obvious from reactants/products, or simply state the transformation clearly.
- Mention specific quantities or conditions if provided in the input text.
- Warn about safety if hazardous materials are obvious (e.g. "Note: Uses KCN, handle with care").

Input Data:
{synthesis_plan}

Output Format:
### Synthesis of [Target Name/SMILES]
[Introduction sentence]

### Step 1: [Transformation Summary]
[Procedure text]
...
### Safety Notes
[General safety warnings]
"""
        prompt = ChatPromptTemplate.from_template(template)
        self.chain = prompt | self.llm | StrOutputParser()

    def get_runnable(self):
        """
        Returns a Runnable that accepts a string (JSON dump of steps) or a dict.
        Supports streaming.
        """
        return RunnableLambda(self._prepare_input) | self.chain

    def _prepare_input(self, input_data):
        # Allow input to be a dict (the full state) or just the string
        if isinstance(input_data, dict):
            # If it's the state dict from run_retrosynthesis
            if "reactions" in input_data:
                 plan_str = json.dumps(input_data, indent=2)
            else:
                 plan_str = str(input_data)
        else:
            plan_str = str(input_data)
        return {"synthesis_plan": plan_str}
