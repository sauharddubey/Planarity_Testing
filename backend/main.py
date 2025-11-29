from fastapi import FastAPI, Body, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from typing import List, Dict, Any
import hashlib
import asyncio
from concurrent.futures import ProcessPoolExecutor
from app.worker import process_graph_task

app = FastAPI()

# Configure CORS
origins = [
    "http://localhost:3000",
    "http://127.0.0.1:3000",
    "http://localhost:5173",
    "http://127.0.0.1:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global cache for graph analysis results
# Key: SHA256 hash of input string
# Value: Result dictionary
RESULT_CACHE: Dict[str, Any] = {}

# Global executor
executor = ProcessPoolExecutor()

@app.on_event("shutdown")
def shutdown_event():
    executor.shutdown()

@app.get("/")
def read_root():
    return {"Hello": "World"}

from fastapi.responses import StreamingResponse
import json

# Define a wrapper to keep track of index and input for process_graph_task
def process_graph_task_with_index(index, inp):
    res = process_graph_task(inp)
    return index, inp, res

@app.post("/process-batch")
async def process_batch(inputs: List[str] = Body(...)):
    """
    Process a batch of graph strings.
    Streams results as NDJSON (Newline Delimited JSON).
    """
    async def result_generator():
        # Use a set to track hashes of inputs already processed or in cache
        processed_hashes = set()
        
        # 1. Cache Hits
        misses = []
        for i, inp in enumerate(inputs):
            h = hashlib.sha256(inp.encode("utf-8")).hexdigest()
            
            if h in RESULT_CACHE:
                if h not in processed_hashes: # Avoid yielding duplicates if same input appears multiple times
                    yield json.dumps({"index": i, "result": RESULT_CACHE[h]}) + "\n"
                    processed_hashes.add(h)
            else:
                misses.append((i, inp, h)) # Store index, input, and its hash for misses
        
        # 2. Process Cache Misses
        if misses:
            loop = asyncio.get_running_loop()
            
            tasks = []
            # Filter out duplicate misses based on hash before creating tasks
            unique_misses_map = {} # hash -> (index, inp)
            for i, inp, h in misses:
                if h not in unique_misses_map:
                    unique_misses_map[h] = (i, inp)
                # If a hash is already in unique_misses_map, it means we've seen this input before.
                # We still need to process it for the original index, but we only need one worker task per unique input.
            
            # Create tasks for unique inputs
            for h, (original_index, unique_inp) in unique_misses_map.items():
                tasks.append(
                    loop.run_in_executor(executor, process_graph_task_with_index, original_index, unique_inp)
                )
            
            # Map hashes back to all original indices that requested them
            hash_to_indices = {}
            for i, inp, h in misses:
                hash_to_indices.setdefault(h, []).append(i)

            for coro in asyncio.as_completed(tasks):
                # The result from process_graph_task_with_index is (original_index, inp, res)
                # We use the 'inp' from the result to re-calculate the hash, ensuring consistency
                # and allowing us to look up all indices that requested this specific input.
                _, inp_from_worker, res = await coro
                
                # Update Cache
                h = hashlib.sha256(inp_from_worker.encode("utf-8")).hexdigest()
                RESULT_CACHE[h] = res
                
                # Yield result for all indices that requested this input
                if h in hash_to_indices:
                    for idx in hash_to_indices[h]:
                        yield json.dumps({"index": idx, "result": res}) + "\n"

    return StreamingResponse(result_generator(), media_type="application/x-ndjson")
