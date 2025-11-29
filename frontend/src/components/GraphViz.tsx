import { useEffect, useRef } from 'react'
import * as d3 from 'd3'

interface Node {
    id: number | string
    x: number
    y: number
    label: string
}

interface Edge {
    source: number | string
    target: number | string
    is_conflict: boolean
}

interface GraphVizProps {
    nodes: Node[]
    edges: Edge[]
    width?: number
    height?: number
}

const GraphViz = ({ nodes, edges, width = 600, height = 400 }: GraphVizProps) => {
    const svgRef = useRef<SVGSVGElement>(null)

    useEffect(() => {
        if (!svgRef.current || nodes.length === 0) return

        const svg = d3.select(svgRef.current)
        svg.selectAll("*").remove() // Clear previous render

        // Create a group for zooming/panning
        const g = svg.append("g")

        // Setup Zoom
        const zoom = d3.zoom<SVGSVGElement, unknown>()
            .scaleExtent([0.1, 10])
            .on("zoom", (event) => {
                g.attr("transform", event.transform)
            })

        svg.call(zoom)

        // Calculate bounds
        const xExtent = d3.extent(nodes, d => d.x) as [number, number]
        const yExtent = d3.extent(nodes, d => d.y) as [number, number]

        const padding = 40

        // Create scales to map node coordinates to SVG dimensions
        // We want to fit the graph within [padding, width-padding] and [padding, height-padding]
        // while maintaining aspect ratio.

        const minX = xExtent[0] || 0
        const maxX = xExtent[1] || 0
        const minY = yExtent[0] || 0
        const maxY = yExtent[1] || 0

        const domainWidth = maxX - minX || 1
        const domainHeight = maxY - minY || 1

        const availableWidth = width - 2 * padding
        const availableHeight = height - 2 * padding

        const scaleFactor = Math.min(availableWidth / domainWidth, availableHeight / domainHeight)

        // Center the graph
        const domainCenterX = (minX + maxX) / 2
        const domainCenterY = (minY + maxY) / 2

        const rangeCenterX = width / 2
        const rangeCenterY = height / 2

        // Transform function: input -> scaled and centered output
        const xScale = (x: number) => rangeCenterX + (x - domainCenterX) * scaleFactor
        const yScale = (y: number) => rangeCenterY + (y - domainCenterY) * scaleFactor

        // Draw Edges
        // Sort edges so conflicts are on top
        const sortedEdges = [...edges].sort((a, b) => (a.is_conflict === b.is_conflict) ? 0 : a.is_conflict ? 1 : -1)

        g.selectAll("line")
            .data(sortedEdges)
            .enter()
            .append("line")
            .attr("x1", d => xScale(nodes.find(n => n.id === d.source)?.x || 0))
            .attr("y1", d => yScale(nodes.find(n => n.id === d.source)?.y || 0))
            .attr("x2", d => xScale(nodes.find(n => n.id === d.target)?.x || 0))
            .attr("y2", d => yScale(nodes.find(n => n.id === d.target)?.y || 0))
            .attr("stroke", d => d.is_conflict ? "#ff4444" : "#999")
            .attr("stroke-width", d => d.is_conflict ? 3 : 1)
            .attr("stroke-opacity", d => d.is_conflict ? 1 : 0.3)
            .attr("class", d => d.is_conflict ? "edge-conflict" : "edge-normal")

        // Draw Nodes
        const nodeGroups = g.selectAll(".node")
            .data(nodes)
            .enter()
            .append("g")
            .attr("class", "node")
            .attr("transform", d => `translate(${xScale(d.x)},${yScale(d.y)})`)

        // Node Circles
        nodeGroups.append("circle")
            .attr("r", 5) // Smaller radius for large graphs
            .attr("fill", "#2b2b2b")
            .attr("stroke", "#646cff")
            .attr("stroke-width", 1.5)
            .on("mouseover", function () { d3.select(this).attr("fill", "#646cff"); })
            .on("mouseout", function () { d3.select(this).attr("fill", "#2b2b2b"); })

        // Node Labels (only if < 50 nodes to avoid clutter, or on hover)
        if (nodes.length < 50) {
            nodeGroups.append("text")
                .text(d => d.label || d.id)
                .attr("dy", ".35em")
                .attr("text-anchor", "middle")
                .attr("fill", "white")
                .style("font-size", "8px")
                .style("pointer-events", "none")
                .style("user-select", "none")
        } else {
            // Add simple title for hover
            nodeGroups.append("title").text(d => d.label || d.id)
        }

    }, [nodes, edges, width, height])

    return (
        <div style={{ width: '100%', height: '100%', border: '1px solid #444', borderRadius: '8px', overflow: 'hidden', background: '#111' }}>
            <svg
                ref={svgRef}
                width="100%"
                height="100%"
                viewBox={`0 0 ${width} ${height}`}
                style={{ cursor: 'grab' }}
            />
        </div>
    )
}

export default GraphViz
