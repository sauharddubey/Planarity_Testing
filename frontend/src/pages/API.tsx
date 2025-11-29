import { useState } from 'react';
import { motion } from 'framer-motion';

const API = () => {
    const [activeSection, setActiveSection] = useState('introduction');

    const scrollToSection = (id: string) => {
        setActiveSection(id);
        const element = document.getElementById(id);
        if (element) {
            element.scrollIntoView({ behavior: 'smooth' });
        }
    };

    return (
        <div className="doc-layout">
            {/* Sidebar */}
            <nav className="doc-sidebar">
                <div className="doc-nav-group">
                    <h3>Getting Started</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'introduction' ? 'active' : ''}`}
                        onClick={() => scrollToSection('introduction')}
                    >
                        Introduction
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Endpoints</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'health-check' ? 'active' : ''}`}
                        onClick={() => scrollToSection('health-check')}
                    >
                        Health Check
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'process-batch' ? 'active' : ''}`}
                        onClick={() => scrollToSection('process-batch')}
                    >
                        Process Batch
                    </a>
                </div>
            </nav>

            {/* Main Content */}
            <main className="doc-content">
                <motion.div
                    initial={{ opacity: 0, y: 10 }}
                    animate={{ opacity: 1, y: 0 }}
                    transition={{ duration: 0.4 }}
                >
                    {/* Introduction */}
                    <section id="introduction" className="api-header-section">
                        <h1 className="api-title">API Reference</h1>
                        <p className="api-description">
                            Welcome to the Planarity Testing API. This API allows you to programmatically analyze graphs for planarity,
                            detect Kuratowski subgraphs (K5, K3,3), and visualize chemical compounds.
                        </p>
                        <div style={{ marginTop: '20px', display: 'flex', gap: '10px', alignItems: 'center' }}>
                            <span style={{ color: 'var(--text-secondary)', fontSize: '0.9rem' }}>Base URL:</span>
                            <code style={{
                                background: 'rgba(255,255,255,0.1)',
                                padding: '6px 12px',
                                borderRadius: '6px',
                                fontFamily: 'JetBrains Mono',
                                color: 'var(--accent-primary)'
                            }}>
                                http://127.0.0.1:8000
                            </code>
                        </div>
                    </section>

                    {/* Health Check */}
                    <section id="health-check" className="endpoint-card">
                        <div className="endpoint-header">
                            <span className="method-badge method-get">GET</span>
                            <span className="endpoint-path">/</span>
                        </div>
                        <p style={{ marginBottom: '20px', color: 'var(--text-secondary)' }}>
                            Verifies that the API service is running and responsive.
                        </p>

                        <div className="code-preview">
                            <div className="code-preview-header">
                                <span>Response</span>
                                <span>200 OK</span>
                            </div>
                            <pre>{`{
  "Hello": "World"
}`}</pre>
                        </div>
                    </section>

                    {/* Process Batch */}
                    <section id="process-batch" className="endpoint-card">
                        <div className="endpoint-header">
                            <span className="method-badge method-post">POST</span>
                            <span className="endpoint-path">/process-batch</span>
                        </div>
                        <p style={{ marginBottom: '20px', color: 'var(--text-secondary)' }}>
                            Submits a batch of graphs for parallel processing. The server uses a race condition between
                            a cache lookup and a fresh computation to return the fastest result.
                        </p>

                        <h3>Request Body</h3>
                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Parameter</th>
                                    <th>Description</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>
                                        <span className="param-name">body</span>
                                        <span className="param-type">string[]</span>
                                    </td>
                                    <td className="param-desc">
                                        An array of graph representations. Each string can be:
                                        <ul style={{ margin: '8px 0', paddingLeft: '20px', color: 'var(--text-secondary)' }}>
                                            <li>Edge List (e.g., <code>"0 1\n1 2"</code>)</li>
                                            <li>SMILES String (e.g., <code>"C1CCCCC1"</code>)</li>
                                            <li>Adjacency Matrix</li>
                                        </ul>
                                    </td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="code-preview">
                            <div className="code-preview-header">
                                <span>Example Request</span>
                                <span>application/json</span>
                            </div>
                            <pre>{`[
  "0 1 2 3 4 0 2 4 1 3",          // K5
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  // Caffeine
]`}</pre>
                        </div>

                        <h3 style={{ marginTop: '30px' }}>Response</h3>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '15px' }}>
                            Returns a stream of <strong>NDJSON</strong> (Newline Delimited JSON). Each line corresponds to one input graph.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Field</th>
                                    <th>Description</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>
                                        <span className="param-name">index</span>
                                        <span className="param-type">integer</span>
                                    </td>
                                    <td className="param-desc">The index of the graph in the input array.</td>
                                </tr>
                                <tr>
                                    <td>
                                        <span className="param-name">result</span>
                                        <span className="param-type">object</span>
                                    </td>
                                    <td className="param-desc">
                                        Contains <code>status</code> ("success" or "error") and <code>data</code>.
                                        <br />
                                        <code>data</code> includes <code>is_planar</code>, <code>nodes</code>, <code>edges</code>, and Kuratowski counts.
                                    </td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="code-preview">
                            <div className="code-preview-header">
                                <span>Example Response Stream</span>
                                <span>application/x-ndjson</span>
                            </div>
                            <pre>{`{"index": 0, "result": {"status": "success", "data": {"is_planar": false, "k5_count": 1, ...}}}
{"index": 1, "result": {"status": "success", "data": {"is_planar": true, "nodes": [...], ...}}}`}</pre>
                        </div>
                    </section>
                </motion.div>
            </main>
        </div>
    );
};

export default API;
