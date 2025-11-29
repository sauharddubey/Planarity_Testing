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
                        className={`doc-nav-item ${activeSection === 'system-design' ? 'active' : ''}`}
                        onClick={() => scrollToSection('system-design')}
                    >
                        System Design
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'architecture' ? 'active' : ''}`}
                        onClick={() => scrollToSection('architecture')}
                    >
                        Architecture
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'data-flow' ? 'active' : ''}`}
                        onClick={() => scrollToSection('data-flow')}
                    >
                        Data Flow
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Backend Components</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'backend-main' ? 'active' : ''}`}
                        onClick={() => scrollToSection('backend-main')}
                    >
                        Main Server
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'backend-parser' ? 'active' : ''}`}
                        onClick={() => scrollToSection('backend-parser')}
                    >
                        Graph Parser
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'backend-engine' ? 'active' : ''}`}
                        onClick={() => scrollToSection('backend-engine')}
                    >
                        Analysis Engine
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'backend-worker' ? 'active' : ''}`}
                        onClick={() => scrollToSection('backend-worker')}
                    >
                        Worker Process
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'backend-drug' ? 'active' : ''}`}
                        onClick={() => scrollToSection('backend-drug')}
                    >
                        Drug Discovery
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Frontend Components</h3>
                    <a
                        className={`doc-nav-item ${activeSection === 'frontend-app' ? 'active' : ''}`}
                        onClick={() => scrollToSection('frontend-app')}
                    >
                        App Structure
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'frontend-pages' ? 'active' : ''}`}
                        onClick={() => scrollToSection('frontend-pages')}
                    >
                        Pages
                    </a>
                    <a
                        className={`doc-nav-item ${activeSection === 'frontend-components' ? 'active' : ''}`}
                        onClick={() => scrollToSection('frontend-components')}
                    >
                        Components
                    </a>
                </div>
                <div className="doc-nav-group">
                    <h3>Workflows</h3>
                    <a className={`doc-nav-item ${activeSection === 'planarity-workflow' ? 'active' : ''}`} onClick={() => scrollToSection('planarity-workflow')}>Planarity Testing</a>
                    <a className={`doc-nav-item ${activeSection === 'drug-workflow' ? 'active' : ''}`} onClick={() => scrollToSection('drug-workflow')}>Drug Discovery</a>
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
                        className={`doc-nav-item ${activeSection === 'input-formats' ? 'active' : ''}`}
                        onClick={() => scrollToSection('input-formats')}
                    >
                        Input Formats
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
                <motion.div initial={{ opacity: 0, y: 10 }} animate={{ opacity: 1, y: 0 }} transition={{ duration: 0.4 }}>

                    {/* Overview */}
                    <section id="overview" className="api-header-section">
                        <h1 className="api-title">System Documentation</h1>
                        <p className="api-description" style={{ fontSize: '1.1rem', lineHeight: '1.8' }}>
                            A modern platform for graph planarity testing and molecular analysis. Built with performance,
                            accuracy, and beautiful visualizations in mind.
                        </p>
                        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: '15px', marginTop: '30px' }}>
                            <div className="card" style={{ background: 'rgba(99, 102, 241, 0.1)', border: '1px solid rgba(99, 102, 241, 0.3)', padding: '20px', textAlign: 'center' }}>
                                <div style={{ fontSize: '2rem', marginBottom: '10px' }}>⚡</div>
                                <div style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>Linear Time Algorithms</div>
                            </div>
                            <div className="card" style={{ background: 'rgba(16, 185, 129, 0.1)', border: '1px solid rgba(16, 185, 129, 0.3)', padding: '20px', textAlign: 'center' }}>
                                <div style={{ fontSize: '2rem', marginBottom: '10px' }}>🎨</div>
                                <div style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>2D & 3D Visualization</div>
                            </div>
                            <div className="card" style={{ background: 'rgba(245, 158, 11, 0.1)', border: '1px solid rgba(245, 158, 11, 0.3)', padding: '20px', textAlign: 'center' }}>
                                <div style={{ fontSize: '2rem', marginBottom: '10px' }}>🧬</div>
                                <div style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>Drug Discovery AI</div>
                            </div>
                            <div className="card" style={{ background: 'rgba(236, 72, 153, 0.1)', border: '1px solid rgba(236, 72, 153, 0.3)', padding: '20px', textAlign: 'center' }}>
                                <div style={{ fontSize: '2rem', marginBottom: '10px' }}>🚀</div>
                                <div style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>Parallel Processing</div>
                            </div>
                        </div>
                    </section>

                    {/* System Design */}
                    <section id="system-design" className="endpoint-card">
                        <h2 style={{ fontSize: '2rem', marginBottom: '20px', background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))', WebkitBackgroundClip: 'text', WebkitTextFillColor: 'transparent' }}>
                            System Architecture
                        </h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '30px', fontSize: '1.05rem', lineHeight: '1.7' }}>
                            A modern full-stack application with parallel processing, intelligent caching, and real-time streaming.
                        </p>

                        <div style={{ display: 'grid', gridTemplateColumns: '1fr', gap: '20px' }}>
                            {/* Frontend Layer */}
                            <div className="card" style={{ background: 'rgba(99, 102, 241, 0.1)', border: '1px solid rgba(99, 102, 241, 0.3)', padding: '25px' }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '15px' }}>
                                    <div style={{ fontSize: '2rem' }}>🎨</div>
                                    <h3 style={{ fontSize: '1.3rem', margin: 0, color: 'var(--accent-primary)' }}>Frontend Layer</h3>
                                </div>
                                <p style={{ color: 'var(--text-secondary)', fontSize: '0.95rem', marginBottom: '15px' }}>
                                    React 19 + TypeScript application with 2D/3D visualization capabilities
                                </p>
                                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '10px' }}>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>React 19</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>UI Framework</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>D3.js</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>2D Graphs</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>Three.js</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>3D Rendering</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>Vite</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>Build Tool</span>
                                    </div>
                                </div>
                            </div>

                            {/* Backend Layer */}
                            <div className="card" style={{ background: 'rgba(16, 185, 129, 0.1)', border: '1px solid rgba(16, 185, 129, 0.3)', padding: '25px' }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '15px' }}>
                                    <div style={{ fontSize: '2rem' }}>⚙️</div>
                                    <h3 style={{ fontSize: '1.3rem', margin: 0, color: 'var(--accent-primary)' }}>Backend Layer</h3>
                                </div>
                                <p style={{ color: 'var(--text-secondary)', fontSize: '0.95rem', marginBottom: '15px' }}>
                                    FastAPI server with async processing and parallel computation
                                </p>
                                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '10px' }}>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>FastAPI</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>Web Framework</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>NetworkX</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>Graph Algorithms</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>RDKit</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>Chemistry</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>Uvicorn</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>ASGI Server</span>
                                    </div>
                                </div>
                            </div>

                            {/* AI & External */}
                            <div className="card" style={{ background: 'rgba(245, 158, 11, 0.1)', border: '1px solid rgba(245, 158, 11, 0.3)', padding: '25px' }}>
                                <div style={{ display: 'flex', alignItems: 'center', gap: '12px', marginBottom: '15px' }}>
                                    <div style={{ fontSize: '2rem' }}>🤖</div>
                                    <h3 style={{ fontSize: '1.3rem', margin: 0, color: 'var(--accent-primary)' }}>AI & External Services</h3>
                                </div>
                                <p style={{ color: 'var(--text-secondary)', fontSize: '0.95rem', marginBottom: '15px' }}>
                                    Google Gemini 2.0 Flash for intelligent molecular analysis
                                </p>
                                <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(150px, 1fr))', gap: '10px' }}>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>Gemini API</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>AI Responses</span>
                                    </div>
                                    <div style={{ padding: '10px', background: 'rgba(255,255,255,0.05)', borderRadius: '6px', fontSize: '0.85rem' }}>
                                        <strong>Async Calls</strong><br />
                                        <span style={{ color: 'var(--text-tertiary)' }}>Non-blocking</span>
                                    </div>
                                </div>
                            </div>
                        </div>

                        <div className="code-preview" style={{ background: '#1a1a2e', padding: '20px', marginTop: '30px' }}>
                            <pre style={{ fontSize: '0.75rem', lineHeight: '1.4' }}>{`
┌─────────────────────────────────────────────────────────────────────────────┐
│                           SYSTEM ARCHITECTURE                                │
└─────────────────────────────────────────────────────────────────────────────┘

┌──────────────────────────────────────────────────────────────────────────────┐
│                              FRONTEND LAYER                                   │
│  ┌────────────────────────────────────────────────────────────────────────┐  │
│  │  React Application (TypeScript + Vite)                                 │  │
│  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐                │  │
│  │  │   Home.tsx   │  │ DrugDiscovery│  │   API.tsx    │                │  │
│  │  │  (Main UI)   │  │   .tsx       │  │   Docs.tsx   │                │  │
│  │  └──────┬───────┘  └──────┬───────┘  └──────────────┘                │  │
│  │         │                  │                                            │  │
│  │  ┌──────▼──────────────────▼────────────────────────────────────────┐  │  │
│  │  │              Shared Components                                    │  │  │
│  │  │  • FileUploader  • GraphViz  • MoleculeViewer                   │  │  │
│  │  │  • Navbar        • StarBackground                                │  │  │
│  │  └───────────────────────────────────────────────────────────────────┘  │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
└────────────────────────────────────┬─────────────────────────────────────────┘
                                     │ HTTP/REST API
                                     │ (axios, fetch)
┌────────────────────────────────────▼─────────────────────────────────────────┐
│                              BACKEND LAYER                                    │
│  ┌────────────────────────────────────────────────────────────────────────┐  │
│  │  FastAPI Server (Python 3.x + Uvicorn)                                │  │
│  │  ┌──────────────────────────────────────────────────────────────────┐ │  │
│  │  │  main.py - API Endpoints & Request Routing                       │ │  │
│  │  │  • GET  /                    (Health Check)                      │ │  │
│  │  │  • POST /process-batch       (Graph Analysis - Streaming)        │ │  │
│  │  │  • POST /drug-discovery/chat (Molecular Analysis + AI)           │ │  │
│  │  └────────────┬─────────────────────────────────────────────────────┘ │  │
│  │               │                                                        │  │
│  │  ┌────────────▼──────────────────────────────────────────────────────┐ │  │
│  │  │  Middleware & Processing Layer                                    │ │  │
│  │  │  ┌──────────────┐  ┌──────────────┐  ┌──────────────────────┐   │ │  │
│  │  │  │   parser.py  │  │  engine.py   │  │  drug_discovery.py   │   │ │  │
│  │  │  │  (Input      │  │  (Graph      │  │  (Molecular          │   │ │  │
│  │  │  │   Parsing)   │  │   Analysis)  │  │   Analysis)          │   │ │  │
│  │  │  └──────┬───────┘  └──────┬───────┘  └──────┬───────────────┘   │ │  │
│  │  │         │                  │                  │                    │ │  │
│  │  │  ┌──────▼──────────────────▼──────────────────▼────────────────┐ │ │  │
│  │  │  │              worker.py (Process Pool)                        │ │ │  │
│  │  │  │  • Parallel Task Execution                                   │ │ │  │
│  │  │  │  • Cache Race Mechanism                                      │ │ │  │
│  │  │  └──────────────────────────────────────────────────────────────┘ │ │  │
│  │  └───────────────────────────────────────────────────────────────────┘ │  │
│  └────────────────────────────────────────────────────────────────────────┘  │
└────────────────────────────────────┬─────────────────────────────────────────┘
                                     │
                    ┌────────────────┴────────────────┐
                    │                                 │
┌───────────────────▼──────────┐   ┌────────────────▼──────────────┐
│   CORE LIBRARIES             │   │   EXTERNAL SERVICES            │
│  • NetworkX (Graph Algo)     │   │  • Google Gemini API           │
│  • RDKit (Chemistry)         │   │    (AI-Powered Responses)      │
│  • NumPy (Numerical)         │   │                                │
└──────────────────────────────┘   └────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                           DATA STORAGE                                       │
│  • In-Memory Cache (SHA-256 keyed results)                                  │
│  • LocalStorage (Frontend state persistence)                                │
└─────────────────────────────────────────────────────────────────────────────┘
`}</pre>
                        </div>
                    </section>

                    {/* Data Flow */}
                    <section id="data-flow" className="endpoint-card">
                        <h2 style={{ fontSize: '2rem', marginBottom: '20px', background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))', WebkitBackgroundClip: 'text', WebkitTextFillColor: 'transparent' }}>
                            Data Flow
                        </h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '30px', fontSize: '1.05rem', lineHeight: '1.7' }}>
                            Follow the journey of data from user input to beautiful visualization.
                        </p>

                        {/* Graph Analysis Flow */}
                        <h3 style={{ fontSize: '1.3rem', marginBottom: '20px', color: 'var(--accent-primary)' }}>📊 Graph Analysis Pipeline</h3>

                        <div style={{ display: 'flex', flexDirection: 'column', gap: '15px', marginBottom: '40px' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>1</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>User Input</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Edge list, SMILES, JSON, or matrix → Stored in React state
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>2</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>HTTP Request</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        POST /process-batch → Array of graph strings sent to backend
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>3</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>Cache Race</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        SHA-256 hash → Race: Cache lookup vs Fresh computation → Winner returns first
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>4</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>Parsing</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Auto-detect format → Convert to NetworkX Graph → Preserve node labels
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>5</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>Analysis</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        LR Algorithm → Planarity test O(n) → Extract certificate → Generate layout
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>6</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>Streaming Response</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        NDJSON stream → Progressive updates → Frontend renders as results arrive
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-primary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>7</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(16, 185, 129, 0.1)', border: '1px solid rgba(16, 185, 129, 0.3)' }}>
                                    <strong style={{ color: 'var(--accent-primary)' }}>Visualization</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Planar → 2D SVG (D3.js) | Non-planar → 3D WebGL (Three.js) | Interactive controls
                                    </p>
                                </div>
                            </div>
                        </div>

                        {/* Drug Discovery Flow */}
                        <h3 style={{ fontSize: '1.3rem', marginBottom: '20px', color: 'var(--accent-primary)' }}>🧬 Drug Discovery Pipeline</h3>

                        <div style={{ display: 'flex', flexDirection: 'column', gap: '15px' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-secondary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>1</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-secondary)' }}>SMILES Input</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Single molecule or batch array → User query about properties
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-secondary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>2</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-secondary)' }}>RDKit Analysis</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Parse SMILES → Calculate MW, LogP, HBD, HBA, TPSA → Lipinski Rule check
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-secondary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>3</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-secondary)' }}>3D Structure</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        ETKDG algorithm → Add hydrogens → Generate 3D coordinates → CPK coloring
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-secondary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>4</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                                    <strong style={{ color: 'var(--accent-secondary)' }}>AI Response</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Gemini 2.0 Flash → Context: molecular properties → Scientific interpretation
                                    </p>
                                </div>
                            </div>

                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px' }}>
                                <div style={{ minWidth: '40px', height: '40px', borderRadius: '50%', background: 'var(--accent-secondary)', display: 'flex', alignItems: 'center', justifyContent: 'center', fontSize: '1.2rem', fontWeight: 'bold' }}>5</div>
                                <div className="card" style={{ flex: 1, padding: '15px', background: 'rgba(245, 158, 11, 0.1)', border: '1px solid rgba(245, 158, 11, 0.3)' }}>
                                    <strong style={{ color: 'var(--accent-secondary)' }}>Display</strong>
                                    <p style={{ margin: '5px 0 0 0', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                        Chat interface → Property cards → 3D molecule viewer (React Three Fiber)
                                    </p>
                                </div>
                            </div>
                        </div>

                        <div className="code-preview" style={{ background: '#1a1a2e', padding: '20px', marginTop: '40px', display: 'none' }}>
                            <pre style={{ fontSize: '0.7rem', lineHeight: '1.4' }}>{`
┌─────────────────────────────────────────────────────────────────────────────┐
│                        GRAPH ANALYSIS DATA FLOW                              │
└─────────────────────────────────────────────────────────────────────────────┘

USER INPUT
    │
    ├─ Edge List:     "0 1\\n1 2\\n2 0"
    ├─ SMILES:        "C1CCCCC1" (Cyclohexane)
    ├─ JSON:          {"nodes": [...], "edges": [...]}
    ├─ Matrix:        [[0,1],[1,0]]
    └─ File Upload:   .txt, .csv, .json
    │
    ▼
┌───────────────────────────────────────────────────────────────────────────┐
│  FRONTEND: Home.tsx                                                        │
│  • User enters graph data in textarea or uploads file                     │
│  • FileUploader component reads file content                              │
│  • Input stored in React state (inputs[])                                 │
│  • User clicks "Analyze All Graphs"                                       │
└───────────────────────────────┬───────────────────────────────────────────┘
                                │
                                ▼ HTTP POST /process-batch
                                │ Body: ["graph1", "graph2", ...]
┌───────────────────────────────▼───────────────────────────────────────────┐
│  BACKEND: main.py - process_batch()                                        │
│  1. Receives array of graph strings                                       │
│  2. For each input, spawns race_for_input() coroutine                     │
│  3. Streams results as NDJSON (Newline Delimited JSON)                    │
└───────────────────────────────┬───────────────────────────────────────────┘
                                │
                    ┌───────────┴───────────┐
                    │                       │
                    ▼                       ▼
        ┌───────────────────┐   ┌───────────────────┐
        │  CACHE LOOKUP     │   │  COMPUTE TASK     │
        │  (async)          │   │  (ProcessPool)    │
        │                   │   │                   │
        │  SHA-256 hash     │   │  worker.py        │
        │  of input         │   │  process_graph    │
        │  → RESULT_CACHE   │   │  _task()          │
        └─────────┬─────────┘   └─────────┬─────────┘
                  │                       │
                  │  asyncio.wait(FIRST_COMPLETED)
                  │                       │
                  └───────────┬───────────┘
                              │ Winner returns first
                              ▼
                    ┌─────────────────────┐
                    │  GRAPH PROCESSING   │
                    └─────────┬───────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  PARSER: parser.py - GraphParser.parse()                                    │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  Input Detection Logic (Sequential Fallthrough):                      │ │
│  │  1. JSON?          → nx.node_link_graph()                             │ │
│  │  2. Matrix?        → nx.from_numpy_array()                            │ │
│  │  3. SMILES?        → RDKit: Chem.MolFromSmiles()                      │ │
│  │     • Atoms → Nodes (with 'label' attribute)                          │ │
│  │     • Bonds → Edges                                                   │ │
│  │  4. Edge List?     → nx.parse_edgelist()                              │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│  Output: NetworkX Graph object                                             │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  ENGINE: engine.py - analyze_graph()                                        │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  STEP 1: Planarity Test                                               │ │
│  │  • nx.check_planarity(G, counterexample=True)                         │ │
│  │  • Returns: (is_planar: bool, certificate: Graph/Embedding)           │ │
│  │                                                                        │ │
│  │  STEP 2: Layout Generation                                            │ │
│  │  IF planar:                                                            │ │
│  │    • nx.planar_layout(G) → Deterministic 2D coordinates               │ │
│  │    • Certificate = PlanarEmbedding (rotation system)                  │ │
│  │  ELSE (non-planar):                                                    │ │
│  │    • nx.spring_layout(G) → Force-directed layout                      │ │
│  │    • Certificate = Kuratowski subgraph (K5 or K3,3)                   │ │
│  │    • Mark conflict edges (edges in Kuratowski subgraph)               │ │
│  │                                                                        │ │
│  │  STEP 3: Serialize Output                                             │ │
│  │  • Nodes: [{id, x, y, label}, ...]                                    │ │
│  │  • Edges: [{source, target, is_conflict}, ...]                        │ │
│  │  • Certificate: Embedding or Kuratowski edges                         │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│  Output: {is_planar, nodes, edges, certificate}                            │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  RESPONSE: Streamed back to frontend                                        │
│  Format: NDJSON (one JSON object per line)                                 │
│  {"index": 0, "result": {"status": "success", "data": {...}}}              │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼ Received by fetch() stream reader
┌─────────────────────────────────────────────────────────────────────────────┐
│  FRONTEND: Home.tsx - processBatch()                                        │
│  • Reads NDJSON stream line by line                                        │
│  • Updates results[] state array at corresponding index                    │
│  • Triggers re-render for each completed graph                             │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  VISUALIZATION: GraphViz.tsx                                                │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  Decision: Planar or Non-Planar?                                      │ │
│  │                                                                        │ │
│  │  IF PLANAR (no conflict edges):                                       │ │
│  │    → 2D Rendering with D3.js                                          │ │
│  │      • SVG-based visualization                                        │ │
│  │      • Force simulation (optional, or use provided coords)            │ │
│  │      • Zoom/Pan controls                                              │ │
│  │      • CPK coloring for atoms (if molecular graph)                    │ │
│  │                                                                        │ │
│  │  IF NON-PLANAR (has conflict edges):                                  │ │
│  │    → 3D Rendering with react-force-graph-3d                           │ │
│  │      • WebGL-based 3D visualization                                   │ │
│  │      • Conflict edges highlighted in red                              │ │
│  │      • Interactive camera controls                                    │ │
│  │      • Animated rotation                                              │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│  Output: Interactive graph visualization                                   │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│                    DRUG DISCOVERY DATA FLOW                                  │
└─────────────────────────────────────────────────────────────────────────────┘

USER INPUT (DrugDiscovery.tsx)
    │
    ├─ SMILES String:  "CC(=O)OC1=CC=CC=C1C(=O)O" (Aspirin)
    ├─ SMILES Array:   ["C1CCCCC1", "c1ccccc1"]
    └─ User Query:     "What are the drug-like properties?"
    │
    ▼ HTTP POST /drug-discovery/chat
    │ Body: {smiles: "...", message: "..."}
┌───────────────────────────────────────────────────────────────────────────┐
│  BACKEND: main.py - drug_discovery_chat()                                  │
│  1. Receives SMILES (single or batch) + user message                      │
│  2. Offloads analysis to ProcessPoolExecutor                              │
│  3. Calls Gemini API for AI response                                      │
└───────────────────────────────┬───────────────────────────────────────────┘
                                │
                    ┌───────────┴───────────┐
                    │                       │
                    ▼                       ▼
        ┌───────────────────┐   ┌───────────────────┐
        │  MOLECULAR        │   │  AI RESPONSE      │
        │  ANALYSIS         │   │  GENERATION       │
        │  (CPU-bound)      │   │  (I/O-bound)      │
        └─────────┬─────────┘   └─────────┬─────────┘
                  │                       │
                  ▼                       │
┌─────────────────────────────────────────────────────────────────────────────┐
│  DRUG DISCOVERY: drug_discovery.py - analyze_molecule()                    │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  STEP 1: Parse SMILES                                                 │ │
│  │  • RDKit: Chem.MolFromSmiles(smiles)                                  │ │
│  │  • Validates molecular structure                                      │ │
│  │                                                                        │ │
│  │  STEP 2: Calculate Properties                                         │ │
│  │  • Molecular Weight (Descriptors.MolWt)                               │ │
│  │  • LogP - Lipophilicity (Descriptors.MolLogP)                         │ │
│  │  • H-Bond Donors (Descriptors.NumHDonors)                             │ │
│  │  • H-Bond Acceptors (Descriptors.NumHAcceptors)                       │ │
│  │  • TPSA - Topological Polar Surface Area                              │ │
│  │  • Rotatable Bonds (Descriptors.NumRotatableBonds)                    │ │
│  │                                                                        │ │
│  │  STEP 3: Lipinski Rule of 5 Check                                     │ │
│  │  • MW ≤ 500 Da                                                         │ │
│  │  • LogP ≤ 5                                                            │ │
│  │  • HBD ≤ 5                                                             │ │
│  │  • HBA ≤ 10                                                            │ │
│  │  • Violations ≤ 1 → Drug-like                                          │ │
│  │                                                                        │ │
│  │  STEP 4: Generate 3D Structure                                        │ │
│  │  • Add hydrogens: Chem.AddHs(mol)                                     │ │
│  │  • Embed 3D coords: AllChem.EmbedMolecule(ETKDG)                      │ │
│  │  • Extract atom positions (x, y, z)                                   │ │
│  │  • Extract bond connections                                           │ │
│  │  • Apply CPK coloring                                                 │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
│  Output: {valid, smiles, properties, lipinski, structure_3d}              │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  AI RESPONSE: drug_discovery.py - generate_llm_response()                  │
│  • Constructs prompt with molecular analysis context                       │
│  • Calls Google Gemini API (gemini-2.0-flash)                              │
│  • Returns scientific interpretation and recommendations                   │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  FRONTEND: DrugDiscovery.tsx                                                │
│  • Displays AI response in chat interface                                  │
│  • Shows molecular properties in analysis card                             │
│  • Renders 3D structure in MoleculeViewer component                        │
└─────────────────────────────────┬───────────────────────────────────────────┘
                                  │
                                  ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  3D VISUALIZATION: MoleculeViewer.tsx                                       │
│  • React Three Fiber (Three.js wrapper)                                    │
│  • Atoms rendered as spheres (CPK colored)                                 │
│  • Bonds rendered as cylinders                                             │
│  • Auto-rotation and orbit controls                                        │
└─────────────────────────────────────────────────────────────────────────────┘
`}</pre>
                        </div>
                    </section>

                    {/* Architecture */}
                    <section id="architecture" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Parallel Cache Race Architecture</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            The system uses a <strong>Race-to-Result</strong> pattern to ensure low latency.
                        </p>

                        <div className="code-preview" style={{ background: '#1a1a2e', padding: '20px' }}>
                            <pre style={{ fontSize: '0.75rem', lineHeight: '1.4' }}>{`
┌─────────────────────────────────────────────────────────────────────────────┐
│                    CACHE RACE MECHANISM (Mermaid-style)                      │
└─────────────────────────────────────────────────────────────────────────────┘

                        Request Arrives
                              │
                              ▼
                    ┌─────────────────┐
                    │  Compute Hash   │
                    │  SHA-256(input) │
                    └────────┬────────┘
                             │
                ┌────────────┴────────────┐
                │                         │
                ▼                         ▼
    ┌───────────────────┐     ┌───────────────────┐
    │   Cache Task      │     │   Compute Task    │
    │   (async)         │     │   (ProcessPool)   │
    │                   │     │                   │
    │  Check in-memory  │     │  1. Parse graph   │
    │  RESULT_CACHE     │     │  2. Analyze       │
    │  dict[hash]       │     │  3. Generate viz  │
    └─────────┬─────────┘     └─────────┬─────────┘
              │                         │
              │   asyncio.wait(         │
              │   FIRST_COMPLETED)      │
              │                         │
              └────────────┬────────────┘
                           │
                    ┌──────▼──────┐
                    │   Winner?   │
                    └──────┬──────┘
                           │
            ┌──────────────┴──────────────┐
            │                             │
            ▼                             ▼
    ┌───────────────┐           ┌───────────────┐
    │  Cache Hit    │           │  Compute Done │
    │  Return       │           │  Update Cache │
    │  immediately  │           │  Return       │
    └───────┬───────┘           └───────┬───────┘
            │                           │
            └───────────┬───────────────┘
                        │
                        ▼
                ┌───────────────┐
                │  Stream to    │
                │  Frontend     │
                │  (NDJSON)     │
                └───────────────┘

Benefits:
• Sub-millisecond response for cached graphs
• No blocking on repeated requests
• Automatic deduplication of identical inputs
• Graceful handling of concurrent requests
`}</pre>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginTop: '20px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>Implementation Details</h3>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)' }}>
                                For every graph submitted, the backend spawns two concurrent asynchronous tasks:
                            </p>
                            <ol style={{ paddingLeft: '20px', color: 'var(--text-secondary)', margin: '10px 0' }}>
                                <li style={{ marginBottom: '8px' }}>
                                    <strong>Cache Lookup:</strong> Checks a SHA-256 hash-keyed in-memory store for previous results. This is nearly instantaneous for local dict access.
                                </li>
                                <li style={{ marginBottom: '8px' }}>
                                    <strong>Compute Task:</strong> Offloads the heavy graph analysis to a <code>ProcessPoolExecutor</code> to avoid blocking the main event loop. This runs in a separate Python process.
                                </li>
                                <li>
                                    <strong>Race Winner:</strong> The system uses <code>asyncio.wait(FIRST_COMPLETED)</code> to return the result from whichever task finishes first, ensuring minimal latency.
                                </li>
                            </ol>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)', marginTop: '10px' }}>
                                If the cache misses, the compute task result is automatically stored for future requests. Identical inputs (same hash) share the same compute task, preventing redundant work.
                            </p>
                        </div>
                    </section>

                    {/* Backend Components */}
                    <section id="backend-main" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Backend: main.py (FastAPI Server)</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            The main entry point for the backend API. Handles HTTP requests, CORS, and coordinates all processing tasks.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Component</th>
                                    <th>Description</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">FastAPI App</span></td>
                                    <td className="param-desc">ASGI web framework with automatic OpenAPI docs, async support, and type validation</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">CORS Middleware</span></td>
                                    <td className="param-desc">Allows cross-origin requests from frontend (localhost:3000, localhost:5173)</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">ProcessPoolExecutor</span></td>
                                    <td className="param-desc">Global executor for CPU-bound tasks. Spawns separate Python processes to avoid GIL contention</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">RESULT_CACHE</span></td>
                                    <td className="param-desc">In-memory dictionary storing SHA-256 hashed results. Key: hash(input), Value: analysis result</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">race_for_input()</span></td>
                                    <td className="param-desc">Async coroutine that races cache lookup vs compute task. Returns fastest result</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">StreamingResponse</span></td>
                                    <td className="param-desc">Streams NDJSON results as they complete, enabling progressive UI updates</td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="code-preview" style={{ marginTop: '20px' }}>
                            <div className="code-preview-header">
                                <span>Key Functions</span>
                            </div>
                            <pre>{`# Endpoints:
GET  /                      → Health check
POST /process-batch         → Graph analysis (streaming)
POST /drug-discovery/chat   → Molecular analysis + AI

# Cache Race Logic:
async def race_for_input(index, inp):
    h = hashlib.sha256(inp.encode()).hexdigest()
    
    # Deduplicate compute tasks
    if h not in unique_computes:
        unique_computes[h] = loop.run_in_executor(
            executor, process_graph_task, inp
        )
    
    cache_task = asyncio.create_task(async_check_cache(h))
    compute_task = asyncio.ensure_future(unique_computes[h])
    
    done, pending = await asyncio.wait(
        [cache_task, compute_task],
        return_when=asyncio.FIRST_COMPLETED
    )
    
    # Return winner's result
    if cache_task in done and cache_task.result():
        return cache_task.result()
    else:
        result = await compute_task
        RESULT_CACHE[h] = result
        return result`}</pre>
                        </div>
                    </section>

                    <section id="backend-parser" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Backend: parser.py (Graph Parser)</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Intelligent input parser that detects and converts multiple graph formats into NetworkX Graph objects.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Format</th>
                                    <th>Detection Method</th>
                                    <th>Parsing Strategy</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">JSON</span></td>
                                    <td className="param-desc">Starts with {'"{"'}</td>
                                    <td className="param-desc">
                                        • Checks for nodes/edges format<br />
                                        • Falls back to nx.node_link_graph()<br />
                                        • Preserves node attributes
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Adjacency Matrix</span></td>
                                    <td className="param-desc">Starts with {'"[["'}</td>
                                    <td className="param-desc">
                                        • Parses as 2D array<br />
                                        • Converts to NumPy array<br />
                                        • Uses nx.from_numpy_array()
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">SMILES</span></td>
                                    <td className="param-desc">No spaces + contains C/N/O</td>
                                    <td className="param-desc">
                                        • RDKit: Chem.MolFromSmiles()<br />
                                        • Atoms → Nodes (with 'label')<br />
                                        • Bonds → Edges
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Edge List</span></td>
                                    <td className="param-desc">Fallback (default)</td>
                                    <td className="param-desc">
                                        • Space/tab separated pairs<br />
                                        • nx.parse_edgelist()<br />
                                        • Supports labeled nodes
                                    </td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="code-preview" style={{ marginTop: '20px' }}>
                            <div className="code-preview-header">
                                <span>Example Inputs</span>
                            </div>
                            <pre>{`# Edge List
"0 1\\n1 2\\n2 0"

# SMILES (Benzene)
"c1ccccc1"

# JSON
{"nodes": [{"id": 0}, {"id": 1}], "edges": [{"source": 0, "target": 1}]}

# Adjacency Matrix
"[[0, 1, 0], [1, 0, 1], [0, 1, 0]]"`}</pre>
                        </div>
                    </section>

                    <section id="backend-engine" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Backend: engine.py (Analysis Engine)</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Core graph analysis engine using NetworkX algorithms. Performs planarity testing and generates visualization layouts.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Function</th>
                                    <th>Algorithm</th>
                                    <th>Output</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">check_planarity()</span></td>
                                    <td className="param-desc">Boyer-Myrvold O(n)</td>
                                    <td className="param-desc">
                                        • is_planar: boolean<br />
                                        • certificate: PlanarEmbedding or Kuratowski subgraph
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">planar_layout()</span></td>
                                    <td className="param-desc">Combinatorial embedding</td>
                                    <td className="param-desc">
                                        • Deterministic 2D coordinates<br />
                                        • No edge crossings<br />
                                        • Used for planar graphs
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">spring_layout()</span></td>
                                    <td className="param-desc">Fruchterman-Reingold</td>
                                    <td className="param-desc">
                                        • Force-directed layout<br />
                                        • Minimizes edge crossings<br />
                                        • Used for non-planar graphs
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Conflict Detection</span></td>
                                    <td className="param-desc">Kuratowski subgraph extraction</td>
                                    <td className="param-desc">
                                        • Identifies K5 or K3,3 subdivisions<br />
                                        • Marks edges as is_conflict: true<br />
                                        • Provides non-planarity certificate
                                    </td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="code-preview" style={{ marginTop: '20px' }}>
                            <div className="code-preview-header">
                                <span>Output Format</span>
                            </div>
                            <pre>{`{
  "is_planar": false,
  "nodes": [
    {"id": 0, "x": 0.5, "y": 0.3, "label": "C"},
    {"id": 1, "x": 0.2, "y": 0.8, "label": "N"}
  ],
  "edges": [
    {"source": 0, "target": 1, "is_conflict": false},
    {"source": 1, "target": 2, "is_conflict": true}
  ],
  "certificate": {
    "type": "Kuratowski Subgraph",
    "edges": [[0, 1], [1, 2], ...]
  }
}`}</pre>
                        </div>
                    </section>

                    <section id="backend-worker" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Backend: worker.py (Process Worker)</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Worker function executed in separate processes via ProcessPoolExecutor. Coordinates parsing and analysis.
                        </p>

                        <div className="code-preview">
                            <div className="code-preview-header">
                                <span>Worker Flow</span>
                            </div>
                            <pre>{`def process_graph_task(input_str: str):
    try:
        # 1. Parse input
        graph = GraphParser.parse(input_str)
        
        # 2. Analyze graph
        result = analyze_graph(graph)
        
        # 3. Return success
        return {"status": "success", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# Executed in separate Python process
# Avoids GIL contention
# Enables true parallelism for CPU-bound tasks`}</pre>
                        </div>
                    </section>

                    <section id="backend-drug" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Backend: drug_discovery.py (Molecular Analysis)</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Comprehensive molecular analysis using RDKit and AI-powered insights via Google Gemini.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Property</th>
                                    <th>Calculation Method</th>
                                    <th>Significance</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">Molecular Weight</span></td>
                                    <td className="param-desc">Descriptors.MolWt()</td>
                                    <td className="param-desc">Drug-like: ≤ 500 Da (Lipinski)</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">LogP</span></td>
                                    <td className="param-desc">Descriptors.MolLogP()</td>
                                    <td className="param-desc">Lipophilicity, membrane permeability (≤ 5)</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">H-Bond Donors</span></td>
                                    <td className="param-desc">Descriptors.NumHDonors()</td>
                                    <td className="param-desc">Affects solubility (≤ 5)</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">H-Bond Acceptors</span></td>
                                    <td className="param-desc">Descriptors.NumHAcceptors()</td>
                                    <td className="param-desc">Affects solubility (≤ 10)</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">TPSA</span></td>
                                    <td className="param-desc">Descriptors.TPSA()</td>
                                    <td className="param-desc">Polar surface area, BBB permeability</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Rotatable Bonds</span></td>
                                    <td className="param-desc">Descriptors.NumRotatableBonds()</td>
                                    <td className="param-desc">Molecular flexibility</td>
                                </tr>
                            </tbody>
                        </table>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginTop: '20px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>3D Structure Generation</h3>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)' }}>
                                Uses RDKit's ETKDG (Experimental-Torsion Knowledge Distance Geometry) algorithm:
                            </p>
                            <ol style={{ paddingLeft: '20px', color: 'var(--text-secondary)', margin: '10px 0', fontSize: '0.9rem' }}>
                                <li>Add explicit hydrogens: <code>Chem.AddHs(mol)</code></li>
                                <li>Generate 3D coordinates: <code>AllChem.EmbedMolecule(mol, AllChem.ETKDG())</code></li>
                                <li>Extract atom positions from conformer</li>
                                <li>Apply CPK coloring (C=gray, O=red, N=blue, etc.)</li>
                                <li>Extract bond connections and types</li>
                            </ol>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginTop: '20px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>AI Integration (Google Gemini)</h3>
                            <p style={{ fontSize: '0.95rem', lineHeight: '1.6', color: 'var(--text-secondary)' }}>
                                Model: <code>gemini-2.0-flash</code>
                            </p>
                            <p style={{ fontSize: '0.9rem', lineHeight: '1.6', color: 'var(--text-secondary)', marginTop: '10px' }}>
                                System Prompt: Acts as expert pharmaceutical scientist, interprets molecular properties,
                                suggests optimizations, compares batch results, and provides drug-likeness assessments.
                            </p>
                        </div>
                    </section>

                    {/* Frontend Components */}
                    <section id="frontend-app" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Frontend: App Structure</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            React 19 application built with TypeScript, Vite, and modern UI libraries.
                        </p>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Technology</th>
                                    <th>Purpose</th>
                                    <th>Version</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">React</span></td>
                                    <td className="param-desc">UI framework with hooks and concurrent features</td>
                                    <td className="param-desc">19.2.0</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">TypeScript</span></td>
                                    <td className="param-desc">Type-safe JavaScript with compile-time checks</td>
                                    <td className="param-desc">5.9.3</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Vite</span></td>
                                    <td className="param-desc">Fast build tool with HMR (Hot Module Replacement)</td>
                                    <td className="param-desc">7.2.4</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">React Router</span></td>
                                    <td className="param-desc">Client-side routing (/, /api, /drug-discovery, etc.)</td>
                                    <td className="param-desc">7.9.6</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Framer Motion</span></td>
                                    <td className="param-desc">Animation library for smooth transitions</td>
                                    <td className="param-desc">12.23.24</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">D3.js</span></td>
                                    <td className="param-desc">Data visualization for 2D graphs</td>
                                    <td className="param-desc">7.9.0</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Three.js</span></td>
                                    <td className="param-desc">WebGL 3D rendering engine</td>
                                    <td className="param-desc">0.181.2</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">React Three Fiber</span></td>
                                    <td className="param-desc">React renderer for Three.js</td>
                                    <td className="param-desc">9.4.0</td>
                                </tr>
                            </tbody>
                        </table>
                    </section>

                    <section id="frontend-pages" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Frontend: Pages</h2>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>Home.tsx</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                Main graph analysis interface with multi-tab support, file upload, and real-time streaming results.
                            </p>
                            <ul style={{ paddingLeft: '20px', color: 'var(--text-secondary)', fontSize: '0.85rem', marginTop: '10px' }}>
                                <li>Hero section with scroll-to-dashboard</li>
                                <li>Multi-graph tabs with status indicators (✅ planar, ❌ non-planar, ⚠️ error)</li>
                                <li>Preset graphs (K5, K3,3, Caffeine)</li>
                                <li>FileUploader integration for batch processing</li>
                                <li>NDJSON streaming with progressive updates</li>
                                <li>LocalStorage persistence for inputs and results</li>
                                <li>Certificate download functionality</li>
                            </ul>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>DrugDiscovery.tsx</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                AI-powered molecular analysis with chat interface and 3D visualization.
                            </p>
                            <ul style={{ paddingLeft: '20px', color: 'var(--text-secondary)', fontSize: '0.85rem', marginTop: '10px' }}>
                                <li>Chat-based interaction with AI scientist</li>
                                <li>SMILES input (single or batch JSON array)</li>
                                <li>Real-time molecular property display</li>
                                <li>3D molecule viewer with rotation</li>
                                <li>Batch analysis comparison</li>
                                <li>Stress test button (100 concurrent requests)</li>
                                <li>Lipinski Rule of 5 validation</li>
                            </ul>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>API.tsx & Docs.tsx</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                Comprehensive API reference and system documentation with interactive navigation.
                            </p>
                        </div>
                    </section>

                    <section id="frontend-components" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Frontend: Reusable Components</h2>

                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Component</th>
                                    <th>Technology</th>
                                    <th>Features</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">FileUploader</span></td>
                                    <td className="param-desc">React + FileReader API</td>
                                    <td className="param-desc">
                                        • Drag & drop support<br />
                                        • Multiple file upload<br />
                                        • Accepts .txt, .csv, .json<br />
                                        • Framer Motion animations
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">GraphViz</span></td>
                                    <td className="param-desc">D3.js + react-force-graph-3d</td>
                                    <td className="param-desc">
                                        • Auto-detects planar vs non-planar<br />
                                        • 2D: SVG with zoom/pan (D3)<br />
                                        • 3D: WebGL rendering (Three.js)<br />
                                        • CPK atom coloring<br />
                                        • Conflict edge highlighting<br />
                                        • Force simulation
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">MoleculeViewer</span></td>
                                    <td className="param-desc">React Three Fiber</td>
                                    <td className="param-desc">
                                        • 3D molecular structure<br />
                                        • Atoms as spheres (CPK colored)<br />
                                        • Bonds as cylinders<br />
                                        • Auto-rotation<br />
                                        • Orbit controls<br />
                                        • Environment lighting
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Navbar</span></td>
                                    <td className="param-desc">React Router</td>
                                    <td className="param-desc">
                                        • Active route highlighting<br />
                                        • Responsive design<br />
                                        • Brand logo
                                    </td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">StarBackground</span></td>
                                    <td className="param-desc">Three.js + maath</td>
                                    <td className="param-desc">
                                        • Animated star field<br />
                                        • 5000 particles<br />
                                        • Slow rotation<br />
                                        • Fixed background layer
                                    </td>
                                </tr>
                            </tbody>
                        </table>
                    </section>

                    {/* Planarity Algorithms */}
                    <section id="boyer-myrvold" className="endpoint-card">
                        <h2 style={{ fontSize: '2rem', marginBottom: '20px', background: 'linear-gradient(135deg, var(--accent-primary), var(--accent-secondary))', WebkitBackgroundClip: 'text', WebkitTextFillColor: 'transparent' }}>
                            Planarity Testing Algorithms
                        </h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '30px', fontSize: '1.05rem', lineHeight: '1.7' }}>
                            Three major algorithms have shaped the field of planarity testing. Each represents a breakthrough in computational efficiency and practical implementation.
                        </p>

                        {/* Hopcroft-Tarjan */}
                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '25px', padding: '25px' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px', marginBottom: '15px' }}>
                                <div style={{ fontSize: '2.5rem' }}>🏛️</div>
                                <div>
                                    <h3 style={{ fontSize: '1.4rem', margin: 0, color: 'var(--accent-primary)' }}>Hopcroft-Tarjan (HT) Algorithm</h3>
                                    <div style={{ fontSize: '0.85rem', color: 'var(--text-tertiary)', marginTop: '5px' }}>1974 • The Pioneer</div>
                                </div>
                            </div>

                            <p style={{ color: 'var(--text-secondary)', lineHeight: '1.7', marginBottom: '15px' }}>
                                The first linear-time planarity testing algorithm. Uses depth-first search (DFS) to build a palm tree,
                                then processes vertices in reverse DFS order to construct a planar embedding.
                            </p>

                            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '15px', marginTop: '20px' }}>
                                <div style={{ padding: '15px', background: 'rgba(16, 185, 129, 0.1)', borderRadius: '8px', border: '1px solid rgba(16, 185, 129, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>✅ The Good</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>First O(n) algorithm</li>
                                        <li>Theoretically elegant</li>
                                        <li>Produces planar embedding</li>
                                        <li>Well-studied and proven</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(245, 158, 11, 0.1)', borderRadius: '8px', border: '1px solid rgba(245, 158, 11, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>⚠️ The Bad</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Complex implementation</li>
                                        <li>Difficult to understand</li>
                                        <li>Large constant factors</li>
                                        <li>Hard to debug</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(239, 68, 68, 0.1)', borderRadius: '8px', border: '1px solid rgba(239, 68, 68, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>❌ The Worst</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Rarely used in practice</li>
                                        <li>No Kuratowski subgraph extraction</li>
                                        <li>Slower than newer methods</li>
                                        <li>Maintenance nightmare</li>
                                    </ul>
                                </div>
                            </div>

                            <div style={{ marginTop: '15px', padding: '12px', background: 'rgba(99, 102, 241, 0.1)', borderRadius: '6px', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                <strong>Time Complexity:</strong> O(n) where n = vertices • <strong>Space:</strong> O(n)
                            </div>
                        </div>

                        {/* Boyer-Myrvold */}
                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '25px', padding: '25px' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px', marginBottom: '15px' }}>
                                <div style={{ fontSize: '2.5rem' }}>⚡</div>
                                <div>
                                    <h3 style={{ fontSize: '1.4rem', margin: 0, color: 'var(--accent-primary)' }}>Boyer-Myrvold (BM) Algorithm</h3>
                                    <div style={{ fontSize: '0.85rem', color: 'var(--text-tertiary)', marginTop: '5px' }}>2004 • The Practical Choice</div>
                                </div>
                            </div>

                            <p style={{ color: 'var(--text-secondary)', lineHeight: '1.7', marginBottom: '15px' }}>
                                A simplified linear-time algorithm that uses edge addition instead of vertex addition.
                                Processes edges incrementally while maintaining planarity, making it much easier to implement and understand.
                            </p>

                            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '15px', marginTop: '20px' }}>
                                <div style={{ padding: '15px', background: 'rgba(16, 185, 129, 0.1)', borderRadius: '8px', border: '1px solid rgba(16, 185, 129, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>✅ The Good</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Still O(n) linear time</li>
                                        <li>Much simpler than HT</li>
                                        <li>Easier to implement</li>
                                        <li>Extracts Kuratowski subgraphs</li>
                                        <li>Better constant factors</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(245, 158, 11, 0.1)', borderRadius: '8px', border: '1px solid rgba(245, 158, 11, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>⚠️ The Bad</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Still complex for beginners</li>
                                        <li>Requires careful bookkeeping</li>
                                        <li>Not as fast as LR in practice</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(239, 68, 68, 0.1)', borderRadius: '8px', border: '1px solid rgba(239, 68, 68, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>❌ The Worst</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Superseded by LR algorithm</li>
                                        <li>More memory overhead than LR</li>
                                        <li>Edge-centric approach less intuitive</li>
                                    </ul>
                                </div>
                            </div>

                            <div style={{ marginTop: '15px', padding: '12px', background: 'rgba(99, 102, 241, 0.1)', borderRadius: '6px', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                <strong>Time Complexity:</strong> O(n) where n = vertices • <strong>Space:</strong> O(n) • <strong>Used by:</strong> Some older NetworkX versions
                            </div>
                        </div>

                        {/* Left-Right (LR) */}
                        <div className="card" style={{ background: 'linear-gradient(135deg, rgba(16, 185, 129, 0.15), rgba(99, 102, 241, 0.15))', border: '2px solid var(--accent-primary)', marginBottom: '25px', padding: '25px' }}>
                            <div style={{ display: 'flex', alignItems: 'center', gap: '15px', marginBottom: '15px' }}>
                                <div style={{ fontSize: '2.5rem' }}>🚀</div>
                                <div>
                                    <h3 style={{ fontSize: '1.4rem', margin: 0, color: 'var(--accent-primary)' }}>Left-Right (LR) Planarity Algorithm</h3>
                                    <div style={{ fontSize: '0.85rem', color: 'var(--text-tertiary)', marginTop: '5px' }}>2008 • The Modern Standard ⭐</div>
                                </div>
                            </div>

                            <div style={{ padding: '12px', background: 'rgba(16, 185, 129, 0.2)', borderRadius: '6px', marginBottom: '15px', border: '1px solid rgba(16, 185, 129, 0.4)' }}>
                                <strong style={{ color: 'var(--accent-primary)' }}>🎯 This is what NetworkX uses!</strong>
                            </div>

                            <p style={{ color: 'var(--text-secondary)', lineHeight: '1.7', marginBottom: '15px' }}>
                                The state-of-the-art planarity testing algorithm by Brandes & others. Uses a left-right partition
                                approach with conflict pairs, making it the fastest and most elegant solution in practice.
                            </p>

                            <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '15px', marginTop: '20px' }}>
                                <div style={{ padding: '15px', background: 'rgba(16, 185, 129, 0.15)', borderRadius: '8px', border: '1px solid rgba(16, 185, 129, 0.4)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>✅ The Good</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>O(n) with best constants</li>
                                        <li>Simplest implementation</li>
                                        <li>Fastest in practice</li>
                                        <li>Clean, readable code</li>
                                        <li>Minimal memory usage</li>
                                        <li>Produces PlanarEmbedding</li>
                                        <li>Extracts Kuratowski subgraphs</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(245, 158, 11, 0.1)', borderRadius: '8px', border: '1px solid rgba(245, 158, 11, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>⚠️ The Bad</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Requires understanding DFS</li>
                                        <li>Conflict pair concept takes time to grasp</li>
                                        <li>Less literature than HT</li>
                                    </ul>
                                </div>
                                <div style={{ padding: '15px', background: 'rgba(239, 68, 68, 0.1)', borderRadius: '8px', border: '1px solid rgba(239, 68, 68, 0.3)' }}>
                                    <div style={{ fontSize: '1.2rem', marginBottom: '8px' }}>❌ The Worst</div>
                                    <ul style={{ fontSize: '0.85rem', color: 'var(--text-secondary)', paddingLeft: '20px', margin: 0 }}>
                                        <li>Honestly? Not much!</li>
                                        <li>It's the best we have</li>
                                        <li>Industry standard for a reason</li>
                                    </ul>
                                </div>
                            </div>

                            <div style={{ marginTop: '15px', padding: '12px', background: 'rgba(99, 102, 241, 0.15)', borderRadius: '6px', fontSize: '0.9rem', color: 'var(--text-secondary)' }}>
                                <strong>Time Complexity:</strong> O(n) where n = vertices • <strong>Space:</strong> O(n) • <strong>Used by:</strong> NetworkX (current), OGDF, Boost Graph Library
                            </div>
                        </div>

                        {/* Comparison Table */}
                        <div style={{ marginTop: '30px' }}>
                            <h3 style={{ fontSize: '1.3rem', marginBottom: '15px', color: 'var(--accent-primary)' }}>Quick Comparison</h3>
                            <table className="param-table">
                                <thead>
                                    <tr>
                                        <th>Algorithm</th>
                                        <th>Year</th>
                                        <th>Complexity</th>
                                        <th>Implementation</th>
                                        <th>Speed (Practice)</th>
                                        <th>Status</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <tr>
                                        <td><span className="param-name">Hopcroft-Tarjan</span></td>
                                        <td className="param-desc">1974</td>
                                        <td className="param-desc">O(n)</td>
                                        <td className="param-desc">Very Hard</td>
                                        <td className="param-desc">Slow</td>
                                        <td className="param-desc">Legacy</td>
                                    </tr>
                                    <tr>
                                        <td><span className="param-name">Boyer-Myrvold</span></td>
                                        <td className="param-desc">2004</td>
                                        <td className="param-desc">O(n)</td>
                                        <td className="param-desc">Moderate</td>
                                        <td className="param-desc">Good</td>
                                        <td className="param-desc">Superseded</td>
                                    </tr>
                                    <tr style={{ background: 'rgba(16, 185, 129, 0.1)' }}>
                                        <td><span className="param-name">Left-Right (LR)</span></td>
                                        <td className="param-desc">2008</td>
                                        <td className="param-desc">O(n)</td>
                                        <td className="param-desc">Easy</td>
                                        <td className="param-desc">Fastest ⚡</td>
                                        <td className="param-desc"><strong>Current Standard ⭐</strong></td>
                                    </tr>
                                </tbody>
                            </table>
                        </div>

                        {/* Code Example */}
                        <div style={{ marginTop: '30px' }}>
                            <h3 style={{ fontSize: '1.3rem', marginBottom: '15px', color: 'var(--accent-primary)' }}>How We Use It</h3>
                            <div className="code-preview">
                                <div className="code-preview-header">
                                    <span>backend/app/engine.py</span>
                                    <span>NetworkX (LR Algorithm)</span>
                                </div>
                                <pre>{`import networkx as nx

def analyze_graph(G: nx.Graph) -> dict:
    # NetworkX uses the LR algorithm internally
    is_planar, certificate = nx.check_planarity(G, counterexample=True)
    
    if is_planar:
        # certificate is a PlanarEmbedding (rotation system)
        pos = nx.planar_layout(G)  # O(n) layout from embedding
    else:
        # certificate is a Kuratowski subgraph (K5 or K3,3)
        pos = nx.spring_layout(G)  # Force-directed for visualization
    
    return {
        "is_planar": is_planar,
        "nodes": [...],
        "edges": [...],
        "certificate": certificate
    }`}</pre>
                            </div>
                        </div>
                    </section>

                    {/* Sequence Diagram */}
                    <section id="sequence" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Request Sequence Diagram</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            Visualizing the parallel execution flow for a graph analysis request.
                        </p>

                        <div className="sequence-diagram" style={{ position: 'relative', padding: '40px 20px', background: '#1e1e1e', borderRadius: '8px', overflowX: 'auto' }}>
                            {/* Actors */}
                            <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '20px', minWidth: '600px' }}>
                                <div style={{ width: '100px', textAlign: 'center', fontWeight: 'bold' }}>User</div>
                                <div style={{ width: '100px', textAlign: 'center', fontWeight: 'bold' }}>Backend</div>
                                <div style={{ width: '100px', textAlign: 'center', fontWeight: 'bold' }}>Cache</div>
                                <div style={{ width: '100px', textAlign: 'center', fontWeight: 'bold' }}>Worker</div>
                            </div>

                            {/* Lifelines */}
                            <div style={{ position: 'absolute', top: '70px', bottom: '20px', left: '70px', width: '2px', background: '#444' }}></div>
                            <div style={{ position: 'absolute', top: '70px', bottom: '20px', left: 'calc(25% + 35px)', width: '2px', background: '#444' }}></div>
                            <div style={{ position: 'absolute', top: '70px', bottom: '20px', left: 'calc(50% + 35px)', width: '2px', background: '#444' }}></div>
                            <div style={{ position: 'absolute', top: '70px', bottom: '20px', left: 'calc(75% + 35px)', width: '2px', background: '#444' }}></div>

                            {/* Messages */}
                            <div style={{ position: 'relative', minWidth: '600px', height: '300px' }}>
                                {/* 1. Request */}
                                <div style={{ position: 'absolute', top: '0', left: '70px', width: '25%', borderTop: '2px solid var(--accent-primary)' }}></div>
                                <div style={{ position: 'absolute', top: '-10px', left: '12%', fontSize: '0.8em' }}>POST /process-batch</div>

                                {/* 2. Fork */}
                                <div style={{ position: 'absolute', top: '40px', left: 'calc(25% + 35px)', width: '25%', borderTop: '2px dashed #aaa' }}></div>
                                <div style={{ position: 'absolute', top: '30px', left: '37%', fontSize: '0.8em' }}>Check Cache (Async)</div>

                                <div style={{ position: 'absolute', top: '60px', left: 'calc(25% + 35px)', width: '50%', borderTop: '2px solid #aaa' }}></div>
                                <div style={{ position: 'absolute', top: '50px', left: '50%', fontSize: '0.8em' }}>Submit Task (ProcessPool)</div>

                                {/* 3. Race */}
                                <div style={{ position: 'absolute', top: '100px', left: 'calc(25% + 35px)', right: '0', textAlign: 'center', color: '#e74c3c', fontSize: '0.9em', fontWeight: 'bold' }}>
                                    RACE: First Completed Wins
                                </div>

                                {/* 4. Return */}
                                <div style={{ position: 'absolute', top: '150px', left: 'calc(50% + 35px)', width: '25%', borderTop: '2px dashed #2ecc71' }}></div>
                                <div style={{ position: 'absolute', top: '140px', left: '37%', fontSize: '0.8em', color: '#2ecc71' }}>Result Found (Hit)</div>

                                <div style={{ position: 'absolute', top: '200px', left: '70px', width: '25%', borderTop: '2px solid var(--accent-primary)' }}></div>
                                <div style={{ position: 'absolute', top: '190px', left: '12%', fontSize: '0.8em' }}>Return JSON</div>
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

                    {/* Input Formats */}
                    <section id="input-formats" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Supported Input Formats</h2>
                        <p style={{ color: 'var(--text-secondary)', marginBottom: '20px' }}>
                            The system accepts multiple graph representation formats, automatically detecting and parsing each type.
                        </p>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>1. Edge List (Text)</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '10px' }}>
                                Space or tab-separated node pairs, one edge per line. Supports both numeric and string node IDs.
                            </p>
                            <div className="code-preview">
                                <pre>{`# Numeric nodes
0 1
1 2
2 0

# Named nodes
A B
B C
C A

# Mixed
node1 node2
node2 node3`}</pre>
                            </div>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>2. SMILES (Chemical)</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '10px' }}>
                                Simplified Molecular Input Line Entry System for chemical structures.
                            </p>
                            <div className="code-preview">
                                <pre>{`# Benzene
c1ccccc1

# Aspirin
CC(=O)OC1=CC=CC=C1C(=O)O

# Caffeine
CN1C=NC2=C1C(=O)N(C(=O)N2C)C

# Ethanol
CCO`}</pre>
                            </div>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>3. JSON (Structured)</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '10px' }}>
                                Node-link format with optional attributes.
                            </p>
                            <div className="code-preview">
                                <pre>{`{
  "nodes": [
    {"id": 0, "label": "A"},
    {"id": 1, "label": "B"}
  ],
  "edges": [
    {"source": 0, "target": 1}
  ]
}`}</pre>
                            </div>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)', marginBottom: '15px' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>4. Adjacency Matrix</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '10px' }}>
                                2D array where matrix[i][j] = 1 indicates an edge between nodes i and j.
                            </p>
                            <div className="code-preview">
                                <pre>{`[[0, 1, 1],
 [1, 0, 1],
 [1, 1, 0]]`}</pre>
                            </div>
                        </div>

                        <div className="card" style={{ background: 'rgba(255,255,255,0.02)', border: '1px solid var(--border-color)' }}>
                            <h3 style={{ fontSize: '1.1rem', marginBottom: '10px', color: 'var(--accent-primary)' }}>5. File Upload</h3>
                            <p style={{ fontSize: '0.9rem', color: 'var(--text-secondary)', marginBottom: '10px' }}>
                                Drag & drop or click to upload .txt, .csv, or .json files containing any of the above formats.
                            </p>
                            <ul style={{ paddingLeft: '20px', color: 'var(--text-secondary)', fontSize: '0.85rem' }}>
                                <li>Multiple files supported</li>
                                <li>Batch processing for JSON arrays</li>
                                <li>Automatic format detection per file</li>
                            </ul>
                        </div>
                    </section>

                    {/* Tech Stack */}
                    <section id="tech-stack" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Complete Technology Stack</h2>

                        <h3 style={{ fontSize: '1.3rem', marginBottom: '15px', color: 'var(--accent-primary)' }}>Backend Dependencies</h3>
                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Library</th>
                                    <th>Purpose</th>
                                    <th>Key Features</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">FastAPI</span></td>
                                    <td className="param-desc">Web framework</td>
                                    <td className="param-desc">Async support, automatic OpenAPI docs, type validation</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Uvicorn</span></td>
                                    <td className="param-desc">ASGI server</td>
                                    <td className="param-desc">High-performance async server, WebSocket support</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">NetworkX</span></td>
                                    <td className="param-desc">Graph algorithms</td>
                                    <td className="param-desc">Boyer-Myrvold planarity, layout algorithms, graph operations</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">RDKit</span></td>
                                    <td className="param-desc">Cheminformatics</td>
                                    <td className="param-desc">SMILES parsing, molecular descriptors, 3D structure generation</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">NumPy</span></td>
                                    <td className="param-desc">Numerical computing</td>
                                    <td className="param-desc">Matrix operations, array processing</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Pydantic</span></td>
                                    <td className="param-desc">Data validation</td>
                                    <td className="param-desc">Type checking, request/response models</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Google Generative AI</span></td>
                                    <td className="param-desc">AI integration</td>
                                    <td className="param-desc">Gemini 2.0 Flash model, async API calls</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">python-dotenv</span></td>
                                    <td className="param-desc">Environment config</td>
                                    <td className="param-desc">Secure API key management</td>
                                </tr>
                            </tbody>
                        </table>

                        <h3 style={{ fontSize: '1.3rem', marginTop: '30px', marginBottom: '15px', color: 'var(--accent-primary)' }}>Frontend Dependencies</h3>
                        <table className="param-table">
                            <thead>
                                <tr>
                                    <th>Library</th>
                                    <th>Purpose</th>
                                    <th>Key Features</th>
                                </tr>
                            </thead>
                            <tbody>
                                <tr>
                                    <td><span className="param-name">React 19</span></td>
                                    <td className="param-desc">UI framework</td>
                                    <td className="param-desc">Concurrent features, hooks, virtual DOM</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">TypeScript</span></td>
                                    <td className="param-desc">Type safety</td>
                                    <td className="param-desc">Compile-time checks, IntelliSense, refactoring</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Vite</span></td>
                                    <td className="param-desc">Build tool</td>
                                    <td className="param-desc">Fast HMR, ES modules, optimized builds</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">React Router</span></td>
                                    <td className="param-desc">Routing</td>
                                    <td className="param-desc">Client-side navigation, nested routes</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">D3.js</span></td>
                                    <td className="param-desc">2D visualization</td>
                                    <td className="param-desc">Force simulation, SVG manipulation, zoom/pan</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Three.js</span></td>
                                    <td className="param-desc">3D rendering</td>
                                    <td className="param-desc">WebGL, geometries, materials, lighting</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">React Three Fiber</span></td>
                                    <td className="param-desc">React + Three.js</td>
                                    <td className="param-desc">Declarative 3D, hooks, automatic cleanup</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">React Three Drei</span></td>
                                    <td className="param-desc">3D helpers</td>
                                    <td className="param-desc">OrbitControls, Environment, PerspectiveCamera</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">react-force-graph-3d</span></td>
                                    <td className="param-desc">3D graph viz</td>
                                    <td className="param-desc">Force-directed 3D graphs, WebGL rendering</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Framer Motion</span></td>
                                    <td className="param-desc">Animations</td>
                                    <td className="param-desc">Spring physics, gestures, layout animations</td>
                                </tr>
                                <tr>
                                    <td><span className="param-name">Axios</span></td>
                                    <td className="param-desc">HTTP client</td>
                                    <td className="param-desc">Promise-based, interceptors, request cancellation</td>
                                </tr>
                            </tbody>
                        </table>
                    </section>

                    {/* File Formats */}
                    <section id="formats" className="endpoint-card">
                        <h2 style={{ fontSize: '1.8rem', marginBottom: '15px' }}>Acceptable File Formats</h2>

                        <h3 style={{ fontSize: '1.1rem', marginTop: '20px', color: 'var(--accent-primary)' }}>1. Graph Upload (JSON)</h3>
                        <div className="code-preview">
                            <pre>{`[
  "0 1\\n1 2\\n2 0",           // Edge List String
  {
    "edges": [[0, 1], [1, 2]] // Object with edges array
  }
]`}</pre>
                        </div>

                        <h3 style={{ fontSize: '1.1rem', marginTop: '20px', color: 'var(--accent-primary)' }}>2. Graph Upload (TXT/CSV)</h3>
                        <div className="code-preview">
                            <pre>{`0 1
1 2
2 0
3 4`}</pre>
                        </div>

                        <h3 style={{ fontSize: '1.1rem', marginTop: '20px', color: 'var(--accent-primary)' }}>3. Drug Discovery (JSON)</h3>
                        <div className="code-preview">
                            <pre>{`[
  "CC(=O)OC1=CC=CC=C1C(=O)O",  // Aspirin
  "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" // Caffeine
]`}</pre>
                        </div>
                    </section>

                </motion.div>
            </main>
        </div >
    );
};

export default Docs;
