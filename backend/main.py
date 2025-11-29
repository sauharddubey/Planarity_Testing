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

@app.post("/process-batch")
async def process_batch(inputs: List[str] = Body(...)):
    """
    Process a batch of graph strings.
    Checks cache first, then processes misses in parallel.
    """
    results = [None] * len(inputs)
    misses_indices = []
    misses_inputs = []

    # 1. Check Cache
    for i, inp in enumerate(inputs):
        # Hash the input
        h = hashlib.sha256(inp.encode("utf-8")).hexdigest()
        
        if h in RESULT_CACHE:
            results[i] = RESULT_CACHE[h]
        else:
            misses_indices.append(i)
            misses_inputs.append(inp)

    # 2. Process Misses in Parallel
    if misses_inputs:
        loop = asyncio.get_running_loop()
        tasks = []
        for inp in misses_inputs:
            # run_in_executor runs the task in the process pool
            tasks.append(loop.run_in_executor(executor, process_graph_task, inp))
        
        # Wait for all tasks to complete
        computed_results = await asyncio.gather(*tasks)
        
        # 3. Update Cache and Results
        for i, res in zip(misses_indices, computed_results):
            results[i] = res
            # Update cache
            h = hashlib.sha256(inputs[i].encode("utf-8")).hexdigest()
            RESULT_CACHE[h] = res

    return results
