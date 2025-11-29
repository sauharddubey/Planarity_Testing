const Docs = () => {
    return (
        <div className="page-container">
            <h1 className="page-title">Documentation</h1>

            <div className="grid-2">
                <div className="card">
                    <h3>Boyer-Myrvold Algorithm</h3>
                    <p>
                        Our engine uses the state-of-the-art Boyer-Myrvold algorithm for linear time planarity testing.
                        It operates in O(n) time complexity.
                    </p>
                </div>

                <div className="card">
                    <h3>Kuratowski Subgraphs</h3>
                    <p>
                        When a graph is non-planar, we extract a Kuratowski subgraph (K5 or K3,3 subdivision) as a certificate of non-planarity.
                    </p>
                </div>

                <div className="card">
                    <h3>Parallel Processing</h3>
                    <p>
                        The backend leverages Python's multiprocessing to handle large batches of graphs efficiently, utilizing all available CPU cores.
                    </p>
                </div>

                <div className="card">
                    <h3>Streaming Architecture</h3>
                    <p>
                        Results are streamed in real-time using NDJSON, ensuring the UI remains responsive even during heavy computational loads.
                    </p>
                </div>
            </div>
        </div>
    )
}

export default Docs
