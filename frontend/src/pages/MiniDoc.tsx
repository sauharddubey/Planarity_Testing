import { motion, useInView } from 'framer-motion';
import { useRef } from 'react';
import './MiniDoc.css';

const MiniDoc = () => {
    return (
        <div className="mini-doc-container">
            {/* Hero Section */}
            <HeroSection />
            
            {/* Problem & Solution */}
            <ProblemSolution />
            
            {/* Core Innovation: Parallelization */}
            <ParallelizationSection />
            
            {/* Cache Race Animation */}
            <CacheRaceSection />
            
            {/* Algorithm Deep Dive */}
            <AlgorithmSection />
            
            {/* Technology Stack */}
            <TechnologySection />
            
            {/* Market Potential */}
            <MarketSection />
            
            {/* Conclusion */}
            <ConclusionSection />
        </div>
    );
};

const HeroSection = () => {
    return (
        <section className="mini-hero">
            <motion.div
                initial={{ opacity: 0, y: 50 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.8 }}
                className="mini-hero-content"
            >
                <motion.div
                    className="mini-hero-icon"
                    animate={{ rotate: 360 }}
                    transition={{ duration: 20, repeat: Infinity, ease: "linear" }}
                >
                    🕸️
                </motion.div>
                <h1 className="mini-hero-title">
                    Graph Planarity Testing
                    <span className="mini-hero-subtitle">Reimagined for Scale</span>
                </h1>
                <p className="mini-hero-description">
                    A modern platform combining linear-time algorithms, parallel processing, 
                    and AI-powered molecular analysis
                </p>
                <div className="mini-hero-stats">
                    <div className="stat-item">
                        <div className="stat-value">O(n)</div>
                        <div className="stat-label">Time Complexity</div>
                    </div>
                    <div className="stat-item">
                        <div className="stat-value">100+</div>
                        <div className="stat-label">Concurrent Graphs</div>
                    </div>
                    <div className="stat-item">
                        <div className="stat-value">&lt;1ms</div>
                        <div className="stat-label">Cache Response</div>
                    </div>
                </div>
            </motion.div>
        </section>
    );
};

import { Canvas } from '@react-three/fiber';
import { OrbitControls, Text, Sphere, Line } from '@react-three/drei';
import * as THREE from 'three';

const Graph3D = () => {
    const groupRef = useRef<THREE.Group>(null);

    // Node positions in 3D space
    const problemNodes = [
        { label: 'Sequential', position: [-4, 4, 0], color: '#fca5a5' },
        { label: 'No Cache', position: [-4, 1.5, 0], color: '#fca5a5' },
        { label: 'Redundant', position: [-4, -1.5, 0], color: '#fca5a5' },
        { label: 'Limited', position: [-4, -4, 0], color: '#fca5a5' },
    ];

    const solutionNodes = [
        { label: 'Deduplicated', position: [4, 4, 0], color: '#86efac' },
        { label: '100+ Concurrent', position: [4, 1.5, 0], color: '#86efac' },
        { label: 'Parallel', position: [4, -1.5, 0], color: '#86efac' },
        { label: 'Cache Race', position: [4, -4, 0], color: '#86efac' },
    ];

    // Connections (crossing edges)
    const edges = [
        { from: problemNodes[0].position, to: solutionNodes[2].position }, // Sequential -> Parallel
        { from: problemNodes[1].position, to: solutionNodes[3].position }, // No Cache -> Cache Race
        { from: problemNodes[2].position, to: solutionNodes[0].position }, // Redundant -> Deduplicated
        { from: problemNodes[3].position, to: solutionNodes[1].position }, // Limited -> 100+ Concurrent
    ];

    return (
        <group ref={groupRef}>
            {/* Problem Nodes */}
            {problemNodes.map((node, idx) => (
                <group key={`problem-${idx}`} position={node.position as [number, number, number]}>
                    <Sphere args={[1.0, 32, 32]}>
                        <meshStandardMaterial color={node.color} roughness={0.3} metalness={0.2} />
                    </Sphere>
                    <Text
                        position={[0, 0, 1.05]}
                        fontSize={0.2}
                        color="#ffffff"
                        anchorX="center"
                        anchorY="middle"
                        maxWidth={1.6}
                        textAlign="center"
                    >
                        {node.label}
                    </Text>
                </group>
            ))}

            {/* Solution Nodes */}
            {solutionNodes.map((node, idx) => (
                <group key={`solution-${idx}`} position={node.position as [number, number, number]}>
                    <Sphere args={[1.0, 32, 32]}>
                        <meshStandardMaterial color={node.color} roughness={0.3} metalness={0.2} />
                    </Sphere>
                    <Text
                        position={[0, 0, 1.05]}
                        fontSize={0.2}
                        color="#ffffff"
                        anchorX="center"
                        anchorY="middle"
                        maxWidth={1.6}
                        textAlign="center"
                    >
                        {node.label}
                    </Text>
                </group>
            ))}

            {/* Edges */}
            {edges.map((edge, idx) => (
                <Line
                    key={`edge-${idx}`}
                    points={[edge.from as [number, number, number], edge.to as [number, number, number]]}
                    color="#9ca3af"
                    lineWidth={4}
                />
            ))}
        </group>
    );
};

const ProblemSolution = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">The Gap We Fill</h2>
                <p className="section-description">
                    Transforming traditional limitations into modern capabilities through intelligent architecture
                </p>

                <motion.div
                    className="graph-3d-container"
                    initial={{ scale: 0.8, opacity: 0 }}
                    animate={isInView ? { scale: 1, opacity: 1 } : {}}
                    transition={{ duration: 0.8, delay: 0.2 }}
                >
                    <div className="graph-3d-wrapper">
                        <Canvas camera={{ position: [0, 0, 8], fov: 50 }}>
                            <ambientLight intensity={0.5} />
                            <pointLight position={[10, 10, 10]} intensity={1} />
                            <pointLight position={[-10, -10, -10]} intensity={0.5} />
                            <Graph3D />
                            <OrbitControls
                                enableZoom={true}
                                enablePan={false}
                                autoRotate={false}
                                minDistance={5}
                                maxDistance={15}
                            />
                        </Canvas>
                    </div>
                    <div className="graph-legend">
                        <div className="legend-item">
                            <div className="legend-dot" style={{ background: '#ef4444' }}></div>
                            <span>Traditional Systems</span>
                        </div>
                        <div className="legend-item">
                            <div className="legend-dot" style={{ background: '#10b981' }}></div>
                            <span>Our Platform</span>
                        </div>
                        <div className="legend-note">
                            Drag to rotate • Scroll to zoom
                        </div>
                    </div>
                </motion.div>

                {/* Key Differentiators */}
                <motion.div
                    className="differentiators"
                    initial={{ y: 50, opacity: 0 }}
                    animate={isInView ? { y: 0, opacity: 1 } : {}}
                    transition={{ duration: 0.6, delay: 0.8 }}
                >
                    <h3>Key Technical Innovations</h3>
                    <div className="diff-grid">
                        <div className="diff-item">
                            <div className="diff-number">01</div>
                            <div className="diff-content">
                                <h4>ProcessPoolExecutor</h4>
                                <p>True parallelism bypassing Python's GIL</p>
                            </div>
                        </div>
                        <div className="diff-item">
                            <div className="diff-number">02</div>
                            <div className="diff-content">
                                <h4>Cache Race Mechanism</h4>
                                <p>Async tasks compete for fastest response</p>
                            </div>
                        </div>
                        <div className="diff-item">
                            <div className="diff-number">03</div>
                            <div className="diff-content">
                                <h4>NDJSON Streaming</h4>
                                <p>Progressive results as they complete</p>
                            </div>
                        </div>
                    </div>
                </motion.div>
            </motion.div>
        </section>
    );
};

const ParallelizationSection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section dark-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">Core Innovation: Parallelization</h2>
                <p className="section-description">
                    Process multiple graphs simultaneously using Python's ProcessPoolExecutor
                </p>

                <div className="parallel-visualization">
                    {[1, 2, 3, 4, 5].map((i) => (
                        <motion.div
                            key={i}
                            className="parallel-task"
                            initial={{ x: -50, opacity: 0 }}
                            animate={isInView ? { x: 0, opacity: 1 } : {}}
                            transition={{
                                duration: 0.5,
                                delay: i * 0.1,
                            }}
                        >
                            <div className="task-label">Graph {i}</div>
                            <div className="task-progress-container">
                                <motion.div
                                    className="task-progress"
                                    initial={{ width: 0 }}
                                    animate={isInView ? { width: "100%" } : {}}
                                    transition={{
                                        duration: 1.5,
                                        delay: i * 0.15 + 0.3,
                                        ease: "easeOut"
                                    }}
                                />
                            </div>
                            <div className="task-status">Complete</div>
                        </motion.div>
                    ))}
                </div>

                <div className="feature-grid">
                    <motion.div
                        className="feature-item"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.8 }}
                    >
                        <h4>True Parallelism</h4>
                        <p>Separate Python processes avoid GIL contention</p>
                    </motion.div>
                    <motion.div
                        className="feature-item"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 1.0 }}
                    >
                        <h4>Deduplication</h4>
                        <p>Identical inputs share the same compute task</p>
                    </motion.div>
                    <motion.div
                        className="feature-item"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 1.2 }}
                    >
                        <h4>Streaming Results</h4>
                        <p>NDJSON stream for progressive UI updates</p>
                    </motion.div>
                </div>
            </motion.div>
        </section>
    );
};

const CacheRaceSection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">Cache Race Mechanism</h2>
                <p className="section-description">
                    Two concurrent tasks compete: O(1) cache lookup vs O(n) computation. Fastest wins.
                </p>

                <div className="technical-diagram">
                    {/* Request Input */}
                    <motion.div
                        className="diagram-node input-node"
                        initial={{ scale: 0 }}
                        animate={isInView ? { scale: 1 } : {}}
                        transition={{ duration: 0.5, delay: 0.2 }}
                    >
                        <div className="node-title">Request</div>
                        <div className="node-detail">Graph Input</div>
                    </motion.div>

                    {/* Hash */}
                    <motion.div
                        className="diagram-arrow"
                        initial={{ scaleX: 0 }}
                        animate={isInView ? { scaleX: 1 } : {}}
                        transition={{ duration: 0.4, delay: 0.4 }}
                    >
                        <div className="arrow-label">SHA-256</div>
                    </motion.div>

                    <motion.div
                        className="diagram-node hash-node"
                        initial={{ scale: 0 }}
                        animate={isInView ? { scale: 1 } : {}}
                        transition={{ duration: 0.5, delay: 0.6 }}
                    >
                        <div className="node-title">Hash Key</div>
                        <div className="node-detail">Unique Identifier</div>
                    </motion.div>

                    {/* Split into two paths */}
                    <div className="diagram-split">
                        {/* Cache Path */}
                        <motion.div
                            className="diagram-path cache-path"
                            initial={{ opacity: 0, y: -20 }}
                            animate={isInView ? { opacity: 1, y: 0 } : {}}
                            transition={{ duration: 0.5, delay: 0.8 }}
                        >
                            <div className="path-header">
                                <div className="path-title">Cache Lookup</div>
                                <div className="complexity-badge fast">O(1)</div>
                            </div>
                            <div className="path-steps">
                                <div className="step">1. Check in-memory dict</div>
                                <div className="step">2. Return if exists</div>
                            </div>
                            <div className="path-time">~0.1ms</div>
                        </motion.div>

                        {/* Compute Path */}
                        <motion.div
                            className="diagram-path compute-path"
                            initial={{ opacity: 0, y: 20 }}
                            animate={isInView ? { opacity: 1, y: 0 } : {}}
                            transition={{ duration: 0.5, delay: 1.0 }}
                        >
                            <div className="path-header">
                                <div className="path-title">Fresh Compute</div>
                                <div className="complexity-badge slow">O(n)</div>
                            </div>
                            <div className="path-steps">
                                <div className="step">1. Parse graph structure</div>
                                <div className="step">2. Run LR algorithm</div>
                                <div className="step">3. Generate layout</div>
                                <div className="step">4. Serialize result</div>
                            </div>
                            <div className="path-time">~50-500ms</div>
                        </motion.div>
                    </div>

                    {/* Race Winner */}
                    <motion.div
                        className="diagram-node result-node"
                        initial={{ scale: 0 }}
                        animate={isInView ? { scale: 1 } : {}}
                        transition={{ duration: 0.5, delay: 1.4 }}
                    >
                        <div className="node-title">asyncio.wait(FIRST_COMPLETED)</div>
                        <div className="node-detail">Winner Returns First</div>
                    </motion.div>
                </div>

                {/* Technical Details */}
                <div className="cache-benefits">
                    <div className="benefit-item">
                        <strong>Constant Time Lookup</strong>
                        <p>Dictionary access is O(1) regardless of cache size</p>
                    </div>
                    <div className="benefit-item">
                        <strong>Linear Time Compute</strong>
                        <p>LR planarity algorithm scales with vertex count</p>
                    </div>
                    <div className="benefit-item">
                        <strong>Async Competition</strong>
                        <p>Both tasks run concurrently, fastest wins</p>
                    </div>
                </div>
            </motion.div>
        </section>
    );
};

const AlgorithmSection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section dark-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">Algorithm Evolution</h2>
                <p className="section-description">
                    50 years of research, culminating in the optimal solution
                </p>

                <div className="timeline">
                    <motion.div
                        className="timeline-item"
                        initial={{ x: -100, opacity: 0 }}
                        animate={isInView ? { x: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.2 }}
                    >
                        <div className="timeline-year">1974</div>
                        <div className="timeline-content">
                            <h4>Hopcroft-Tarjan</h4>
                            <p>First O(n) algorithm. Complex, hard to implement.</p>
                            <div className="timeline-badge legacy">Legacy</div>
                        </div>
                    </motion.div>

                    <motion.div
                        className="timeline-item"
                        initial={{ x: -100, opacity: 0 }}
                        animate={isInView ? { x: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.4 }}
                    >
                        <div className="timeline-year">2004</div>
                        <div className="timeline-content">
                            <h4>Boyer-Myrvold</h4>
                            <p>Simplified edge-addition approach. Better, but not best.</p>
                            <div className="timeline-badge superseded">Superseded</div>
                        </div>
                    </motion.div>

                    <motion.div
                        className="timeline-item active"
                        initial={{ x: -100, opacity: 0 }}
                        animate={isInView ? { x: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.6 }}
                    >
                        <div className="timeline-year">2008</div>
                        <div className="timeline-content">
                            <h4>Left-Right (LR)</h4>
                            <p>Optimal O(n) with best constants. Industry standard.</p>
                            <div className="timeline-badge current">Current ⭐</div>
                        </div>
                    </motion.div>
                </div>

                <motion.div
                    className="why-lr-card"
                    initial={{ y: 50, opacity: 0 }}
                    animate={isInView ? { y: 0, opacity: 1 } : {}}
                    transition={{ delay: 0.8 }}
                >
                    <h3>Why Use NetworkX (LR Algorithm)?</h3>
                    <div className="why-lr-grid">
                        <div className="why-item">
                            <div className="why-icon">✅</div>
                            <p><strong>Battle-tested</strong> by thousands of researchers</p>
                        </div>
                        <div className="why-item">
                            <div className="why-icon">⚡</div>
                            <p><strong>Fastest</strong> implementation available</p>
                        </div>
                        <div className="why-item">
                            <div className="why-icon">🔧</div>
                            <p><strong>Maintained</strong> by core Python community</p>
                        </div>
                        <div className="why-item">
                            <div className="why-icon">📚</div>
                            <p><strong>Well-documented</strong> with extensive tests</p>
                        </div>
                        <div className="why-item">
                            <div className="why-icon">🎯</div>
                            <p><strong>Focus on innovation</strong>, not reinventing wheels</p>
                        </div>
                        <div className="why-item">
                            <div className="why-icon">🚀</div>
                            <p><strong>Production-ready</strong> out of the box</p>
                        </div>
                    </div>
                </motion.div>
            </motion.div>
        </section>
    );
};

const TechnologySection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    const technologies = [
        { name: "FastAPI", icon: "⚡", color: "#009688" },
        { name: "React 19", icon: "⚛️", color: "#61DAFB" },
        { name: "NetworkX", icon: "🕸️", color: "#4CAF50" },
        { name: "RDKit", icon: "🧬", color: "#E91E63" },
        { name: "Three.js", icon: "🎨", color: "#000000" },
        { name: "Gemini AI", icon: "🤖", color: "#4285F4" },
    ];

    return (
        <section className="mini-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">Technology Stack</h2>
                <p className="section-description">
                    Modern, performant, and production-ready
                </p>

                <div className="tech-grid">
                    {technologies.map((tech, i) => (
                        <motion.div
                            key={tech.name}
                            className="tech-card"
                            initial={{ scale: 0, rotate: -180 }}
                            animate={isInView ? { scale: 1, rotate: 0 } : {}}
                            transition={{
                                type: "spring",
                                stiffness: 260,
                                damping: 20,
                                delay: i * 0.1,
                            }}
                            whileHover={{ scale: 1.1, rotate: 5 }}
                        >
                            <div className="tech-icon" style={{ fontSize: '3rem' }}>{tech.icon}</div>
                            <div className="tech-name">{tech.name}</div>
                        </motion.div>
                    ))}
                </div>
            </motion.div>
        </section>
    );
};

const MarketSection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section dark-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0 }}
                animate={isInView ? { opacity: 1 } : {}}
                transition={{ duration: 0.6 }}
            >
                <h2 className="mini-section-title">Market Potential</h2>
                <p className="section-description">
                    Powerful API + AI integration = Endless possibilities
                </p>

                <div className="market-grid">
                    <motion.div
                        className="market-card"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.2 }}
                    >
                        <div className="market-icon">🧬</div>
                        <h3>Drug Discovery</h3>
                        <p>AI-powered molecular analysis with Lipinski Rule validation</p>
                        <ul>
                            <li>Pharmaceutical research</li>
                            <li>Chemical compound screening</li>
                            <li>Property prediction</li>
                        </ul>
                    </motion.div>

                    <motion.div
                        className="market-card"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.4 }}
                    >
                        <div className="market-icon">🔬</div>
                        <h3>Research & Academia</h3>
                        <p>Fast, reliable planarity testing for graph theory research</p>
                        <ul>
                            <li>Algorithm benchmarking</li>
                            <li>Educational tools</li>
                            <li>Publication-ready visualizations</li>
                        </ul>
                    </motion.div>

                    <motion.div
                        className="market-card"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.6 }}
                    >
                        <div className="market-icon">🏗️</div>
                        <h3>VLSI & Circuit Design</h3>
                        <p>Verify circuit planarity for PCB layout optimization</p>
                        <ul>
                            <li>Circuit board design</li>
                            <li>Network topology</li>
                            <li>Infrastructure planning</li>
                        </ul>
                    </motion.div>

                    <motion.div
                        className="market-card"
                        initial={{ y: 50, opacity: 0 }}
                        animate={isInView ? { y: 0, opacity: 1 } : {}}
                        transition={{ delay: 0.8 }}
                    >
                        <div className="market-icon">📊</div>
                        <h3>Data Visualization</h3>
                        <p>Automatic layout selection for optimal graph rendering</p>
                        <ul>
                            <li>Network analysis tools</li>
                            <li>Social network visualization</li>
                            <li>Knowledge graphs</li>
                        </ul>
                    </motion.div>
                </div>

                <motion.div
                    className="api-highlight"
                    initial={{ scale: 0.8, opacity: 0 }}
                    animate={isInView ? { scale: 1, opacity: 1 } : {}}
                    transition={{ delay: 1.0 }}
                >
                    <h3>🚀 RESTful API Ready</h3>
                    <p>Easy integration with streaming support and automatic OpenAPI documentation</p>
                </motion.div>
            </motion.div>
        </section>
    );
};

const ConclusionSection = () => {
    const ref = useRef(null);
    const isInView = useInView(ref, { once: true, margin: "-100px" });

    return (
        <section className="mini-section conclusion-section" ref={ref}>
            <motion.div
                initial={{ opacity: 0, y: 50 }}
                animate={isInView ? { opacity: 1, y: 0 } : {}}
                transition={{ duration: 0.8 }}
                className="conclusion-content"
            >
                <motion.div
                    className="conclusion-icon"
                    animate={{ scale: [1, 1.2, 1] }}
                    transition={{ duration: 2, repeat: Infinity }}
                >
                    ✨
                </motion.div>
                <h2>Goals Achieved</h2>
                <div className="achievement-grid">
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>Linear-time O(n) planarity testing</p>
                    </div>
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>Parallel batch processing (100+ graphs)</p>
                    </div>
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>Intelligent cache race mechanism</p>
                    </div>
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>Beautiful 2D/3D visualizations</p>
                    </div>
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>AI-powered drug discovery</p>
                    </div>
                    <div className="achievement-item">
                        <div className="achievement-check">✓</div>
                        <p>Production-ready RESTful API</p>
                    </div>
                </div>

                <motion.div
                    className="final-cta"
                    initial={{ scale: 0 }}
                    animate={isInView ? { scale: 1 } : {}}
                    transition={{ delay: 0.5, type: "spring" }}
                >
                    <h3>Graph Planarity Testing, Perfected.</h3>
                    <p>Fast. Scalable. Beautiful.</p>
                </motion.div>
            </motion.div>
        </section>
    );
};

export default MiniDoc;
