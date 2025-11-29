import { useState } from 'react';
import { motion } from 'framer-motion';

const Docs = () => {
    const [activeSection, setActiveSection] = useState('overview');

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
                    <h3>Guide</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'overview' ? 'active' : ''}`}
                        onClick={() => scrollToSection('overview')}
                    >
                        Overview
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'architecture' ? 'active' : ''}`}
                        onClick={() => scrollToSection('architecture')}
                    >
                        Architecture
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Algorithms</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'boyer-myrvold' ? 'active' : ''}`}
                        onClick={() => scrollToSection('boyer-myrvold')}
                    >
                        Boyer-Myrvold
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'drawing' ? 'active' : ''}`}
                        onClick={() => scrollToSection('drawing')}
                    >
                        Graph Drawing
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'kuratowski' ? 'active' : ''}`}
                        onClick={() => scrollToSection('kuratowski')}
                    >
                        Kuratowski Subgraphs
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Resources</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'tech-stack' ? 'active' : ''}`}
                        onClick={() => scrollToSection('tech-stack')}
                    >
                        Tech Stack
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'references' ? 'active' : ''}`}
                        onClick={() => scrollToSection('references')}
                    >
                        References
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
                    {/* Overview */}
                    <section id="overview" className="api-header-section">
                        <h1 className="api-title">Documentation</h1>
                        <p className="api-description">
                            Deep dive into the algorithms, architecture, and libraries that power the Planarity Testing engine.
                        </p>
                    </section>

                    {/* Architecture */}
                    <section id="architecture" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>System Architecture</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            The application is built for performance and responsiveness, utilizing a parallel processing pipeline that races cache lookups against fresh computations.
                        </p>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>Parallel Cache Race</h3>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)' }}>
                                For every graph submitted, the backend spawns two concurrent asynchronous tasks:
                            </p>
                            <ol style={{ paddingLeft: '20px', color: 'var(--text-secondary)', margin: '10px 0' }}>
                                <li style={{ marginBottom: '8px' }}>
                                    <strong>Cache Lookup:</strong> Checks a SHA-256 hash-keyed in-memory store for previous results.
                                </li>
                                <li>
                                    <strong>Compute Task:</strong> Offloads the heavy graph analysis to a <code>ProcessPoolExecutor</code> to avoid blocking the main event loop.
                                </li>
                            </ol>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)', marginTop: '10px' }}>
                                The system uses <code>asyncio.wait(FIRST_COMPLETED)</code> to return the result from whichever task finishes first, ensuring minimal latency.
                            </p>
                        </div>
                    </section>

                    {/* Boyer-Myrvold */}
                    <section id="boyer-myrvold" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Boyer-Myrvold Algorithm</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            We utilize the <strong>Boyer-Myrvold algorithm</strong> for planarity testing. This is a state-of-the-art algorithm that operates in linear time <code>O(n)</code>, where <code>n</code> is the number of vertices. Unlike the Hopcroft-Tarjan algorithm, Boyer-Myrvold does not require a complex embedding phase and can directly isolate Kuratowski subgraphs.
                        </p>

                        <div className="code-preview">
                            <div className="code-preview-header">
                                <span>Python Implementation</span>
                                <span>networkx</span>
                            </div>
                            <pre>{`import networkx as nx

def check_planarity(G):
    is_planar, certificate = nx.check_planarity(G)
    if is_planar:
        return True, nx.combinatorial_embedding_to_pos(certificate)
    else:
        return False, certificate # Contains Kuratowski subgraph`}</pre>
                        </div>
                    </section>

                    {/* Graph Drawing */}
                    <section id="drawing" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Graph Drawing Algorithms</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Visualizing graphs effectively requires placing nodes in a way that minimizes edge crossings and reveals symmetries. We employ two primary strategies:
                        </p>

                        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '20px' }}>
                            <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>Force-Directed Layout</h3>
                                <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                    Used for <strong>non-planar</strong> graphs and 3D visualizations. Based on the <strong>Fruchterman-Reingold</strong> algorithm, it simulates nodes as charged particles that repel each other and edges as springs that attract connected nodes.
                                </p>
                            </div>
                            <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>Combinatorial Embedding</h3>
                                <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                    Used for <strong>planar</strong> graphs. This deterministic layout ensures that no edges cross by utilizing the planar embedding found by the Boyer-Myrvold algorithm.
                                </p>
                            </div>
                        </div>
                    </section>

                    {/* Kuratowski */}
                    <section id="kuratowski" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Kuratowski Subgraphs</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            According to Kuratowski's theorem, a finite graph is planar if and only if it does not contain a subgraph that is a subdivision of <strong>K<sub>5</sub></strong> (complete graph on 5 vertices) or <strong>K<sub>3,3</sub></strong> (complete bipartite graph on 3+3 vertices).
                        </p>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            When a graph is determined to be non-planar, our engine extracts these specific subgraphs to provide a "certificate" of non-planarity, which is then highlighted in the UI.
                        </p>
                    </section>

                    {/* Tech Stack */}
                    <section id="tech-stack" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Libraries & Tools</h2>
                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Library</th>
                                    <th>Purpose</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td>
                                        <span className="param-name">NetworkX</span>
                                    </td>
                                    <td className="param-desc">
                                        Core graph data structures and algorithms (Boyer-Myrvold, Spring Layout).
                                    </td>
                                </tr>
                                <tr>
                                    <td>
                                        <span className="param-name">RDKit</span>
                                    </td>
                                    <td className="param-desc">
                                        Cheminformatics library used to parse SMILES strings into graph structures for molecular analysis.
                                    </td>
                                </tr>
                                <tr>
                                    <td>
                                        <span className="param-name">FastAPI</span>
                                    </td>
                                    <td className="param-desc">
                                        High-performance async web framework for the backend API and streaming responses.
                                    </td>
                                </tr>
                                <tr>
                                    <td>
                                        <span className="param-name">React Force Graph</span>
                                    </td>
                                    <td className="param-desc">
                                        Frontend library for 2D and 3D force-directed graph rendering using Three.js and D3.
                                    </td>
                                </tr>
                            </tbody>
                        </table>
                    </section>

                    {/* References */}
                    <section id="references" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>References</h2>
                        <ul style={{ listStyle: 'none', padding: 0, color: 'var(--text-secondary)', fontSize: '0.95rem' }}>
                            <li style={{ marginBottom: '15px', paddingLeft: '15px', borderLeft: '2px solid var(--accent-primary)' }}>
                                <strong>[1] J. M. Boyer and W. J. Myrvold</strong>, "On the Cutting Edge: Simplified O(n) Planarity by Edge Addition," <em>Journal of Graph Algorithms and Applications</em>, vol. 8, no. 3, pp. 241-273, 2004.
                            </li>
                            <li style={{ marginBottom: '15px', paddingLeft: '15px', borderLeft: '2px solid var(--accent-primary)' }}>
                                <strong>[2] T. M. J. Fruchterman and E. M. Reingold</strong>, "Graph Drawing by Force-directed Placement," <em>Software: Practice and Experience</em>, vol. 21, no. 11, pp. 1129-1164, 1991.
                            </li>
                            <li style={{ marginBottom: '15px', paddingLeft: '15px', borderLeft: '2px solid var(--accent-primary)' }}>
                                <strong>[3] K. Kuratowski</strong>, "Sur le problème des courbes gauches en topologie," <em>Fundamenta Mathematicae</em>, vol. 15, pp. 271-283, 1930.
                            </li>
                        </ul>
                    </section>

                </motion.div>
            </main>
        </div >
    );
};

export default Docs;
