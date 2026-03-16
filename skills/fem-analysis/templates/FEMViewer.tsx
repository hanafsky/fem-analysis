// FEMViewer.tsx — Interactive 3D FEM result viewer
// Used with: npm create vite@latest viewer -- --template react-ts
// Dependencies: npm install three @react-three/fiber @react-three/drei

import { useMemo, useState, useEffect } from 'react'
import { Canvas } from '@react-three/fiber'
import { OrbitControls, Html } from '@react-three/drei'
import * as THREE from 'three'

interface FEMData {
  vertices: number[][]     // [[x, y], ...]
  triangles: number[][]    // [[i, j, k], ...]
  scalars: {
    name: string
    values: number[]
    min: number
    max: number
  }
}

// Viridis-like colormap
function viridisColor(t: number): [number, number, number] {
  t = Math.max(0, Math.min(1, t))
  const r = Math.max(0, Math.min(1, -0.87 * t * t + 1.52 * t + 0.26))
  const g = Math.max(0, Math.min(1, 0.11 * t * t + 0.65 * t + 0.15))
  const b = Math.max(0, Math.min(1, -1.64 * t * t + 0.92 * t + 0.53))
  return [r, g, b]
}

function FEMMesh({ data }: { data: FEMData }) {
  const geometry = useMemo(() => {
    const geom = new THREE.BufferGeometry()
    const { vertices, triangles, scalars } = data

    const positions = new Float32Array(triangles.length * 3 * 3)
    const colors = new Float32Array(triangles.length * 3 * 3)

    triangles.forEach((tri, i) => {
      tri.forEach((vi, j) => {
        const idx = (i * 3 + j) * 3
        positions[idx] = vertices[vi][0]
        positions[idx + 1] = vertices[vi][1]
        positions[idx + 2] = 0

        const t = (scalars.values[vi] - scalars.min) / (scalars.max - scalars.min || 1)
        const [r, g, b] = viridisColor(t)
        colors[idx] = r
        colors[idx + 1] = g
        colors[idx + 2] = b
      })
    })

    geom.setAttribute('position', new THREE.BufferAttribute(positions, 3))
    geom.setAttribute('color', new THREE.BufferAttribute(colors, 3))
    geom.computeVertexNormals()
    return geom
  }, [data])

  return (
    <mesh geometry={geometry}>
      <meshBasicMaterial vertexColors side={THREE.DoubleSide} />
    </mesh>
  )
}

function ColorBar({ min, max, name }: { min: number; max: number; name: string }) {
  const steps = 5
  return (
    <div style={{
      position: 'absolute', right: 20, top: '50%', transform: 'translateY(-50%)',
      background: 'rgba(255,255,255,0.9)', padding: 10, borderRadius: 8, fontSize: 12
    }}>
      <div style={{ fontWeight: 'bold', marginBottom: 4 }}>{name}</div>
      {Array.from({ length: steps + 1 }, (_, i) => {
        const t = 1 - i / steps
        const val = min + t * (max - min)
        const [r, g, b] = viridisColor(t)
        return (
          <div key={i} style={{ display: 'flex', alignItems: 'center', gap: 4 }}>
            <div style={{
              width: 20, height: 16,
              background: `rgb(${r * 255},${g * 255},${b * 255})`
            }} />
            <span>{val.toExponential(2)}</span>
          </div>
        )
      })}
    </div>
  )
}

export default function FEMViewer({ jsonUrl = '/result.json' }: { jsonUrl?: string }) {
  const [data, setData] = useState<FEMData | null>(null)

  useEffect(() => {
    fetch(jsonUrl)
      .then(r => r.json())
      .then(setData)
      .catch(e => console.error('Failed to load FEM data:', e))
  }, [jsonUrl])

  if (!data) return <div style={{ padding: 20 }}>Loading FEM data...</div>

  // Auto-center and scale
  const xs = data.vertices.map(v => v[0])
  const ys = data.vertices.map(v => v[1])
  const cx = (Math.min(...xs) + Math.max(...xs)) / 2
  const cy = (Math.min(...ys) + Math.max(...ys)) / 2
  const span = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys))

  return (
    <div style={{ width: '100%', height: '100vh', position: 'relative' }}>
      <Canvas camera={{ position: [cx, cy, span * 1.5], fov: 50, near: 0.001, far: 1000 }}>
        <ambientLight intensity={0.8} />
        <FEMMesh data={data} />
        <OrbitControls target={[cx, cy, 0]} />
      </Canvas>
      <ColorBar min={data.scalars.min} max={data.scalars.max} name={data.scalars.name} />
    </div>
  )
}
