import json
import sys
from pathlib import Path
from ord_schema.message_helpers import load_message
from ord_schema.proto import dataset_pb2
from google.protobuf.json_format import MessageToJson

# Set up paths
data_root = Path("../ord-data/data")
output_dir = Path("data")

# Create output directory if it doesn't exist
output_dir.mkdir(exist_ok=True)

# Find all .pb.gz files recursively in the data root
pb_files = list(data_root.rglob("*.pb.gz"))

if not pb_files:
    print(f"No .pb.gz files found in {data_root}")
    sys.exit(0)

print(f"Found {len(pb_files)} collection(s). Starting conversion...")

for input_pb_gz in pb_files:
    # Use the parent directory name or filename as part of the output name
    collection_name = input_pb_gz.stem.split('.')[0] # e.g. 'ord_dataset-xxx'
    # If the parent is a 2-char folder like '00', '0a', we prefer that
    subdir_name = input_pb_gz.parent.name
    output_jsonl = output_dir / f"ord_{subdir_name}_{collection_name}.jsonl"

    print(f"Converting {input_pb_gz}...")
    
    # Load the entire ORD dataset
    try:
        dataset = load_message(str(input_pb_gz), dataset_pb2.Dataset)
    except Exception as e:
        print(f"Failed to load {input_pb_gz}: {e}")
        continue

    # Open output for streaming
    with open(output_jsonl, "w") as out:
        for reaction in dataset.reactions:
            # Convert the reaction protobuf to JSON string
            json_str = MessageToJson(
                reaction,
                including_default_value_fields=False,
                preserving_proto_field_name=True,
            )
            # Parse to Python dict
            reaction_dict = json.loads(json_str)
            # Write each reaction as one line
            out.write(json.dumps(reaction_dict) + "\n")

    print(f"  Saved to {output_jsonl}")

print("\nAll conversions complete.")
