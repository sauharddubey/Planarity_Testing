const API = () => {
    return (
        <div className="page-container">
            <h1 className="page-title">API Reference</h1>
            <div className="card">
                <h2>POST /process-batch</h2>
                <p>Submit a batch of graphs for planarity testing.</p>

                <div className="code-block">
                    <div className="code-header">Request Body</div>
                    <pre>{`[
  "0 1 2 3 4 0 2 4 1 3", // K5
  "0 3 0 4 0 5 1 3 1 4 1 5 2 3 2 4 2 5" // K3,3
]`}</pre>
                </div>

                <div className="code-block">
                    <div className="code-header">Response (NDJSON Stream)</div>
                    <pre>{`{"index": 0, "result": {"is_planar": false, ...}}
{"index": 1, "result": {"is_planar": false, ...}}`}</pre>
                </div>
            </div>
        </div>
    )
}

export default API
