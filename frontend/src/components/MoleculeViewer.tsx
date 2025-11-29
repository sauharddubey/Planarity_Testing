import React, { useRef } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, PerspectiveCamera, Environment } from '@react-three/drei';
import * as THREE from 'three';

interface Atom {
    symbol: string;
    x: number;
    y: number;
    z: number;
    color: string;
}

interface Bond {
    start: { x: number; y: number; z: number };
    end: { x: number; y: number; z: number };
    order: number;
}

interface MoleculeData {
    atoms: Atom[];
    bonds: Bond[];
}

interface MoleculeViewerProps {
    data: MoleculeData;
}

const AtomMesh: React.FC<{ atom: Atom }> = ({ atom }) => {
    // Scale radius based on element (simplified)
    const radius = atom.symbol === 'H' ? 0.25 : 0.4;
    return (
        <mesh position={[atom.x, atom.y, atom.z]}>
            <sphereGeometry args={[radius, 32, 32]} />
            <meshStandardMaterial color={atom.color} roughness={0.3} metalness={0.2} />
        </mesh>
    );
};

const BondMesh: React.FC<{ bond: Bond }> = ({ bond }) => {
    const start = new THREE.Vector3(bond.start.x, bond.start.y, bond.start.z);
    const end = new THREE.Vector3(bond.end.x, bond.end.y, bond.end.z);
    const direction = new THREE.Vector3().subVectors(end, start);
    const length = direction.length();

    // Position is midpoint
    const position = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5);

    // Orientation
    const quaternion = new THREE.Quaternion();
    quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction.clone().normalize());

    // Bond thickness based on order (simplified visualization)
    const radius = 0.1;

    return (
        <mesh position={position} quaternion={quaternion}>
            <cylinderGeometry args={[radius, radius, length, 16]} />
            <meshStandardMaterial color="#CCCCCC" roughness={0.5} metalness={0.1} />
        </mesh>
    );
};

const RotatingMolecule: React.FC<{ data: MoleculeData }> = ({ data }) => {
    const groupRef = useRef<THREE.Group>(null);

    useFrame((_, delta) => {
        if (groupRef.current) {
            groupRef.current.rotation.y += delta * 0.2; // Slow rotation
        }
    });

    // Center the molecule
    const center = new THREE.Vector3();
    if (data.atoms.length > 0) {
        data.atoms.forEach(atom => center.add(new THREE.Vector3(atom.x, atom.y, atom.z)));
        center.divideScalar(data.atoms.length);
    }

    return (
        <group ref={groupRef} position={[-center.x, -center.y, -center.z]}>
            {data.atoms.map((atom, i) => (
                <AtomMesh key={`atom-${i}`} atom={atom} />
            ))}
            {data.bonds.map((bond, i) => (
                <BondMesh key={`bond-${i}`} bond={bond} />
            ))}
        </group>
    );
};

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ data }) => {
    return (
        <div style={{ width: '100%', height: '300px', background: 'rgba(0,0,0,0.2)', borderRadius: '12px', overflow: 'hidden' }}>
            <Canvas>
                <PerspectiveCamera makeDefault position={[0, 0, 10]} />
                <ambientLight intensity={0.5} />
                <pointLight position={[10, 10, 10]} intensity={1} />
                <spotLight position={[-10, -10, -10]} intensity={0.5} />
                <Environment preset="city" />

                <RotatingMolecule data={data} />

                <OrbitControls enableZoom={true} enablePan={true} autoRotate={false} />
            </Canvas>
        </div>
    );
};

export default MoleculeViewer;
