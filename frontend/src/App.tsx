import { useState } from 'react'
import axios from 'axios'
import './App.css'

interface GraphData {
  is_planar: boolean
  nodes: Array<{ id: number | string; x: number; y: number; label: string }>
  edges: Array<{ source: number | string; target: number | string; is_conflict: boolean }>
}

interface Result {
  status: 'success' | 'error'
  data?: GraphData
  message?: string
}

function App() {
  const [input, setInput] = useState<string>('1 2\n2 3\n3 1\n\n---\n\nC-C-O')
  const [results, setResults] = useState<Result[]>([])
  const [loading, setLoading] = useState<boolean>(false)
  const [activeTab, setActiveTab] = useState<number>(0)
  const [error, setError] = useState<string | null>(null)

  const handleAnalyze = async () => {
    setLoading(true)
    setError(null)
    setResults([])
    setActiveTab(0)

    try {
      // Split input by double newlines or "---" delimiter
      // Normalize newlines first
      const normalizedInput = input.replace(/\r\n/g, '\n')

      let inputs: string[] = []
      if (normalizedInput.includes('---')) {
        inputs = normalizedInput.split('---').map(s => s.trim()).filter(s => s.length > 0)
      } else {
        inputs = normalizedInput.split('\n\n').map(s => s.trim()).filter(s => s.length > 0)
      }

      if (inputs.length === 0) {
        setError("Please enter at least one graph.")
        setLoading(false)
        return
      }

      const response = await axios.post('http://localhost:8000/process-batch', inputs)
      setResults(response.data)
    } catch (err) {
      console.error(err)
      setError("Failed to fetch results. Ensure backend is running.")
    } finally {
      setLoading(false)
    }
  }

  return (
    <div className="dashboard-container">
      <div className="left-panel">
        <h2>Input</h2>
        <p style={{ fontSize: '0.9em', color: '#aaa' }}>
          Enter graphs separated by double newlines or "---".<br />
          Formats: Edge List ("1 2"), SMILES ("C-C-O"), Adjacency Matrix ("[[0,1],[1,0]]").
        </p>
        <textarea
          value={input}
          onChange={(e) => setInput(e.target.value)}
          placeholder="Enter graph data here..."
        />
        <button onClick={handleAnalyze} disabled={loading}>
          {loading ? 'Analyzing...' : 'Analyze'}
        </button>
        {error && <div style={{ color: 'red', marginTop: '10px' }}>{error}</div>}
      </div>

      <div className="right-panel">
        <h2>Results</h2>
        {results.length > 0 ? (
          <>
            <div className="tabs">
              {results.map((res, idx) => {
                const isPlanar = res.status === 'success' && res.data?.is_planar
                const isError = res.status === 'error'
                let tabClass = 'tab'
                if (idx === activeTab) tabClass += ' active'
                if (isError) tabClass += ' error'
                else if (isPlanar) tabClass += ' planar'
                else tabClass += ' non-planar'

                return (
                  <div
                    key={idx}
                    className={tabClass}
                    onClick={() => setActiveTab(idx)}
                  >
                    Graph {idx + 1} {isError ? '⚠️' : (isPlanar ? '✅' : '❌')}
                  </div>
                )
              })}
            </div>

            <div className="result-content">
              {results[activeTab] && (
                <div>
                  {results[activeTab].status === 'success' ? (
                    <>
                      <div className={`status-badge ${results[activeTab].data?.is_planar ? 'status-planar' : 'status-non-planar'}`}>
                        {results[activeTab].data?.is_planar ? 'PLANAR' : 'NON-PLANAR'}
                      </div>
                      <div>
                        <strong>Nodes:</strong> {results[activeTab].data?.nodes.length}<br />
                        <strong>Edges:</strong> {results[activeTab].data?.edges.length}
                      </div>
                      <hr style={{ borderColor: '#444', margin: '10px 0' }} />
                      <pre>{JSON.stringify(results[activeTab].data, null, 2)}</pre>
                    </>
                  ) : (
                    <>
                      <div className="status-badge status-error">ERROR</div>
                      <p>{results[activeTab].message}</p>
                    </>
                  )}
                </div>
              )}
            </div>
          </>
        ) : (
          <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', color: '#666' }}>
            No results to display. Run analysis to see results.
          </div>
        )}
      </div>
    </div>
  )
}

export default App
