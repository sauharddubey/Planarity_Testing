import urllib.request
import urllib.error
import json
import time
import concurrent.futures

BASE_URL = "http://127.0.0.1:8000"
ENDPOINT = "/drug-discovery/chat"
ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"

def make_request(i):
    url = f"{BASE_URL}{ENDPOINT}"
    payload = {
        "smiles": ASPIRIN_SMILES,
        "message": f"Request {i}"
    }
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    
    try:
        with urllib.request.urlopen(req) as response:
            if response.status == 200:
                return True, json.loads(response.read().decode())
            return False, None
    except Exception as e:
        return False, str(e)

def test_single_request():
    print("Testing single request...")
    start = time.time()
    success, result = make_request(0)
    duration = time.time() - start
    
    if success:
        print(f"Success! Duration: {duration:.2f}s")
        print(f"Analysis: {json.dumps(result.get('analysis', {}).get('properties'), indent=2)}")
        print(f"Response: {result.get('response')[:100]}...")
        return True
    else:
        print(f"Failed: {result}")
        return False

def test_concurrency(n=20):
    print(f"\nTesting concurrency with {n} requests...")
    start = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=n) as executor:
        results = list(executor.map(make_request, range(n)))
    
    duration = time.time() - start
    success_count = sum(1 for s, _ in results if s)
    print(f"Completed {n} requests in {duration:.2f}s")
    print(f"Success rate: {success_count}/{n}")
    return success_count == n

if __name__ == "__main__":
    if test_single_request():
        test_concurrency(20)
