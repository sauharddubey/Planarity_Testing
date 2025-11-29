import { useEffect, useRef } from 'react'
import * as d3 from 'd3'

interface Node {
    id: number | string
    x: number
    y: number
    label?: string
}

interface Edge {
    source: number | string
    target: number | string
    is_conflict?: boolean
}

interface GraphVizProps {
    nodes: Node[]
    edges: Edge[]
    width?: number
    height?: number
}

const GraphViz = ({ nodes, edges, width = 800, height = 600 }: GraphVizProps) => {
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
            .on("zoom", (event: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
                g.attr("transform", event.transform.toString())
            })

        svg.call(zoom)

        // Calculate bounds
        const xExtent = d3.extent(nodes, (d: Node) => d.x) as [number, number]
        const yExtent = d3.extent(nodes, (d: Node) => d.y) as [number, number]

        const padding = 40

        // Create scales to map node coordinates to SVG dimensions
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
            .attr("x1", (d: Edge) => xScale(nodes.find(n => n.id === d.source)?.x || 0))
            .attr("y1", (d: Edge) => yScale(nodes.find(n => n.id === d.source)?.y || 0))
            .attr("x2", (d: Edge) => xScale(nodes.find(n => n.id === d.target)?.x || 0))
            .attr("y2", (d: Edge) => yScale(nodes.find(n => n.id === d.target)?.y || 0))
            .attr("stroke", (d: Edge) => d.is_conflict ? "#ff4444" : "#999")
            .attr("stroke-width", (d: Edge) => d.is_conflict ? 3 : 1)
            .attr("stroke-opacity", (d: Edge) => d.is_conflict ? 1 : 0.3)
            .attr("class", (d: Edge) => d.is_conflict ? "edge-conflict" : "edge-normal")

        // Draw Nodes
        const nodeGroups = g.selectAll(".node")
            .data(nodes)
            .enter()
            .append("g")
            .attr("class", "node")
            .attr("transform", (d: Node) => `translate(${xScale(d.x)},${yScale(d.y)})`)

        // Node Circles
        nodeGroups.append("circle")
            .attr("r", 6)
            .attr("fill", "#2b2b2b")
            .attr("stroke", "#646cff")
            .attr("stroke-width", 1.5)
            .on("mouseover", function (this: SVGCircleElement) { d3.select(this).attr("fill", "#646cff"); })
            .on("mouseout", function (this: SVGCircleElement) { d3.select(this).attr("fill", "#2b2b2b"); })

        // Node Labels
        // Show labels if they exist (e.g. chemistry atoms) or if graph is small
        nodeGroups.each(function (this: SVGGElement, d: Node) {
            const el = d3.select(this);
            if (d.label) {
                el.append("text")
                    .text(d.label)
                    .attr("dy", ".35em")
                    .attr("text-anchor", "middle")
                    .attr("fill", "white")
                    .style("font-size", "8px")
                    .style("pointer-events", "none")
                    .style("user-select", "none");
            } else if (nodes.length < 50) {
                el.append("text")
                    .text(String(d.id))
                    .attr("dy", ".35em")
                    .attr("text-anchor", "middle")
                    .attr("fill", "white")
                    .style("font-size", "8px")
                    .style("pointer-events", "none")
                    .style("user-select", "none");
            } else {
                el.append("title").text(String(d.id));
            }
        });

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
