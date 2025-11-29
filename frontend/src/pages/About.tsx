const About = () => {
    return (
        <div className="page-container">
            <h1 className="page-title">About PlanarityTest</h1>
            <div className="card">
                <p className="lead">
                    PlanarityTest is a high-performance, modern web application for analyzing graph planarity.
                </p>
                <p>
                    Built with a Python FastAPI backend and a React frontend, it demonstrates the power of modern web technologies combined with efficient algorithms.
                </p>
                <br />
                <p>
                    <strong>Version:</strong> 2.0.0 (Web3 Edition)<br />
                    <strong>Engine:</strong> NetworkX + Custom C++ Extensions (Planned)<br />
                    <strong>License:</strong> MIT
                </p>
            </div>
        </div>
    )
}

export default About
