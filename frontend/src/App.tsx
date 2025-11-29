import { BrowserRouter as Router, Routes, Route } from 'react-router-dom'
import Navbar from './components/Navbar'
import Home from './pages/Home'
import API from './pages/API'
import Docs from './pages/Docs'
import About from './pages/About'
import DrugDiscovery from './pages/DrugDiscovery'
import StarBackground from './components/StarBackground'
import './App.css'

function App() {
  return (
    <Router>
      <div className="app-container">
        <StarBackground />
        <div className="glow-bg" />
        <Navbar />
        <div className="main-content">
          <Routes>
            <Route path="/" element={<Home />} />
            <Route path="/api" element={<API />} />
            <Route path="/docs" element={<Docs />} />
            <Route path="/about" element={<About />} />
            <Route path="/drug-discovery" element={<DrugDiscovery />} />
          </Routes>
        </div>
      </div>
    </Router>
  )
}

export default App
