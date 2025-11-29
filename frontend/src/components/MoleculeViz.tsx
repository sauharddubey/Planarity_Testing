import React, { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';

interface Atom {
    id: number;
    symbol: string;
    x: number;
    y: number;
    color: string;
}

interface Bond {
    source: number | Atom;
    target: number | Atom;
    order: number;
}

interface MoleculeVizProps {
    structure: {
        nodes: Atom[];
        links: Bond[];
    };
    width?: number;
    height?: number;
}

const MoleculeViz: React.FC<MoleculeVizProps> = ({ structure, width = 300, height = 250 }) => {
    const svgRef = useRef<SVGSVGElement>(null);
    const [zoomLevel, setZoomLevel] = useState(1);

    useEffect(() => {
        if (!svgRef.current || !structure) return;

        const svg = d3.select(svgRef.current);
        svg.selectAll("*").remove(); // Clear previous

        const { nodes, links } = structure;

        // RDKit coordinates are centered around 0,0 and usually small (e.g. -5 to 5)
        // We need to scale them to fit the SVG
        const xExtent = d3.extent(nodes, d => d.x) as [number, number];
        const yExtent = d3.extent(nodes, d => d.y) as [number, number];

        // Handle single atom case
        const xMin = xExtent[0] ?? 0;
        const xMax = xExtent[1] ?? 0;
        const yMin = yExtent[0] ?? 0;
        const yMax = yExtent[1] ?? 0;

        const contentWidth = xMax - xMin || 1;
        const contentHeight = yMax - yMin || 1;

        // Add padding
        const padding = 20;

        // Calculate scale to fit
        const scaleX = (width - padding * 2) / contentWidth;
        const scaleY = (height - padding * 2) / contentHeight;
        const scale = Math.min(scaleX, scaleY) * 0.8; // 0.8 factor for extra breathing room

        // Center offset
        const cx = (xMin + xMax) / 2;
        const cy = (yMin + yMax) / 2;

        // Transform function: invert Y because SVG y increases downwards, chemical y increases upwards
        const transformX = (x: number) => (x - cx) * scale + width / 2;
        const transformY = (y: number) => -(y - cy) * scale + height / 2;

        // Group for zoom
        const g = svg.append("g");

        const zoom = d3.zoom<SVGSVGElement, unknown>()
            .scaleExtent([0.5, 5])
            .on("zoom", (event) => {
                g.attr("transform", event.transform.toString());
                setZoomLevel(event.transform.k);
            });

        svg.call(zoom);

        // Render Bonds
        g.append("g")
            .selectAll("line")
            .data(links)
            .join("line")
            .attr("x1", d => transformX(nodes[d.source as number].x))
            .attr("y1", d => transformY(nodes[d.source as number].y))
            .attr("x2", d => transformX(nodes[d.target as number].x))
            .attr("y2", d => transformY(nodes[d.target as number].y))
            .attr("stroke", "#cbd5e1")
            .attr("stroke-width", d => d.order * 2) // Thicker for double bonds
            .attr("opacity", 0.8);

        // Render Atoms
        const atomGroup = g.append("g")
            .selectAll("g")
            .data(nodes)
            .join("g")
            .attr("transform", d => `translate(${transformX(d.x)}, ${transformY(d.y)})`);

        // Atom Circles
        atomGroup.append("circle")
            .attr("r", 12)
            .attr("fill", "#1e293b") // Dark background for contrast
            .attr("stroke", d => d.color)
            .attr("stroke-width", 2);

        // Atom Labels
        atomGroup.append("text")
            .text(d => d.symbol)
            .attr("text-anchor", "middle")
            .attr("dy", "0.35em")
            .attr("fill", d => d.color)
            .attr("font-weight", "bold")
            .attr("font-size", "12px")
            .style("pointer-events", "none");

    }, [structure, width, height]);

    return (
        <div style={{ width, height, background: 'rgba(0,0,0,0.2)', borderRadius: '8px', overflow: 'hidden' }}>
            <svg ref={svgRef} width={width} height={height} style={{ cursor: 'move' }} />
        </div>
    );
};

export default MoleculeViz;
