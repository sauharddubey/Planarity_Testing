import React, { useState, useRef, useEffect } from 'react';
import { motion } from 'framer-motion';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import MoleculeViz from './MoleculeViz';

interface Message {
    role: 'user' | 'assistant';
    content: string;
    analysis?: any; // Result from backend
}

const DrugDiscoveryChat: React.FC = () => {
    const [messages, setMessages] = useState<Message[]>([
        { role: 'assistant', content: "Hello! I'm your AI Drug Discovery Assistant. You can ask me about chemistry or paste SMILES strings (e.g., `CCO`, `c1ccccc1`) for analysis." }
    ]);
    const [input, setInput] = useState('');
    const [isLoading, setIsLoading] = useState(false);
    const messagesEndRef = useRef<HTMLDivElement>(null);

    const scrollToBottom = () => {
        messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
    };

    useEffect(() => {
        scrollToBottom();
    }, [messages]);

    // Heuristic to detect SMILES
    const detectSmiles = (text: string): string[] => {
        // Simple heuristic: 
        // 1. Not a graph edge list (digits and spaces only)
        // 2. Contains typical SMILES chars
        // 3. Length > 1

        const potentialSmiles = text.split(/[\s,]+/).map(t => {
            // Strip surrounding parens/punctuation
            return t.replace(/^[\(\)\.,]+|[\(\)\.,]+$/g, '');
        }).filter(token => {
            if (token.length < 1) return false;
            if (/^\d+$/.test(token)) return false; // Pure numbers
            if (/^[0-9\s]+$/.test(token)) return false; // Graph edge list like "0 1"

            // Check for valid SMILES characters (subset)
            // We allow all letters in the regex initially, but will validate unbracketed ones below
            const validChars = /^[A-Za-z0-9@\.\+\-\[\]\(\)\\\/=#$]+$/;
            if (!validChars.test(token)) return false;

            // Must contain at least one atom symbol
            if (!/[A-Zc-n]/.test(token)) return false;

            // CRITICAL: Filter out tokens containing invalid unbracketed characters.
            // Valid organic subset (plus some common ones): B, C, N, O, P, S, F, Cl, Br, I. 
            // Also aromatic: b, c, n, o, p, s.
            // Everything else must be in brackets (e.g. [Au]).
            // If we see 'd', 'r', 'u', 'g', 'l', 'i', 'k', 'e' outside brackets, it's text.
            // Note: 'l' is in 'Cl', 'r' is in 'Br'. We need to be careful.
            // Strategy: Remove all valid organic tokens and bracketed groups. If anything alphabetic remains, it's invalid.

            let remaining = token.replace(/\[.*?\]/g, ''); // Remove brackets
            // Remove valid 2-letter elements: Cl, Br, In, ... (add more if needed, but for drug discovery organic is key)
            remaining = remaining.replace(/Cl|Br/g, '');
            // Remove valid 1-letter organic elements and aromatics
            remaining = remaining.replace(/[BCNOPSFIbcnops]/g, '');
            // Remove non-alpha chars
            remaining = remaining.replace(/[^a-zA-Z]/g, '');

            if (remaining.length > 0) {
                // Found invalid characters (like 'd', 'u', 'g' from 'drug')
                return false;
            }

            // Filter out common English words (letters only, mixed case or title case)
            if (/^[A-Za-z]+$/.test(token)) {
                const elements = new Set(["C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "B", "c", "n", "o", "p", "s"]);
                if (elements.has(token)) return true;
                if (/^[A-Z]+$/.test(token)) return true; // CCO
                return false;
            }

            // Reject tokens with hyphens unless they look like specific chemical syntax (rare in simple SMILES)
            // "drug-likeness", "3-like" have hyphens.
            if (token.includes('-') && !/\[.*-.*\]/.test(token) && !/[-+=#$]/.test(token)) {
                // If hyphen is not part of charge [-] or bond -, reject.
                // Actually SMILES uses - for single bond, but usually between atoms.
                // Heuristic: if it looks like a word with hyphen
                if (/[a-z]-[a-z]/i.test(token)) return false;
            }

            // Reject "e.g"
            if (token === 'e.g' || token === 'e.g.') return false;

            return true;
        });

        return potentialSmiles;
    };

    const handleSend = async () => {
        if (!input.trim()) return;

        const userMsg = input;
        setInput('');
        setMessages(prev => [...prev, { role: 'user', content: userMsg }]);
        setIsLoading(true);

        try {
            const smiles = detectSmiles(userMsg);

            const payload = {
                message: userMsg,
                smiles: smiles.length > 0 ? smiles : null
            };

            const response = await fetch('http://127.0.0.1:8000/drug-discovery/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });

            const data = await response.json();

            setMessages(prev => [...prev, {
                role: 'assistant',
                content: data.response,
                analysis: data.analysis
            }]);

        } catch (error) {
            console.error("Error sending message:", error);
            setMessages(prev => [...prev, { role: 'assistant', content: "Sorry, I encountered an error processing your request." }]);
        } finally {
            setIsLoading(false);
        }
    };

    const handleKeyDown = (e: React.KeyboardEvent) => {
        if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault();
            handleSend();
        }
    };

    return (
        <div className="chat-container" style={{
            display: 'flex',
            flexDirection: 'column',
            height: '80vh',
            maxWidth: '1000px',
            margin: '0 auto',
            background: 'rgba(30, 41, 59, 0.5)',
            borderRadius: '16px',
            border: '1px solid rgba(255, 255, 255, 0.1)',
            overflow: 'hidden'
        }}>
            {/* Messages Area */}
            <div className="messages-area" style={{
                flex: 1,
                overflowY: 'auto',
                padding: '20px',
                display: 'flex',
                flexDirection: 'column',
                gap: '20px'
            }}>
                {messages.map((msg, idx) => (
                    <motion.div
                        key={idx}
                        initial={{ opacity: 0, y: 10 }}
                        animate={{ opacity: 1, y: 0 }}
                        style={{
                            alignSelf: msg.role === 'user' ? 'flex-end' : 'flex-start',
                            maxWidth: '85%'
                        }}
                    >
                        <div style={{
                            padding: '12px 16px',
                            borderRadius: '12px',
                            background: msg.role === 'user' ? 'var(--accent-primary)' : 'rgba(51, 65, 85, 0.8)',
                            color: '#fff',
                            boxShadow: '0 4px 6px rgba(0,0,0,0.1)'
                        }}>
                            <ReactMarkdown
                                remarkPlugins={[remarkGfm]}
                                components={{
                                    strong: ({ node, ...props }) => <span style={{ fontWeight: 'bold', color: '#fbbf24' }} {...props} />,
                                    em: ({ node, ...props }) => <span style={{ fontStyle: 'italic', color: '#94a3b8' }} {...props} />,
                                    h1: ({ node, ...props }) => <h1 style={{ fontSize: '1.5em', fontWeight: 'bold', margin: '10px 0' }} {...props} />,
                                    h2: ({ node, ...props }) => <h2 style={{ fontSize: '1.3em', fontWeight: 'bold', margin: '8px 0' }} {...props} />,
                                    h3: ({ node, ...props }) => <h3 style={{ fontSize: '1.1em', fontWeight: 'bold', margin: '6px 0' }} {...props} />,
                                    ul: ({ node, ...props }) => <ul style={{ paddingLeft: '20px', margin: '10px 0' }} {...props} />,
                                    li: ({ node, ...props }) => <li style={{ marginBottom: '4px' }} {...props} />,
                                    p: ({ node, ...props }) => <p style={{ margin: '8px 0', lineHeight: '1.5' }} {...props} />,
                                }}
                            >
                                {msg.content}
                            </ReactMarkdown>
                        </div>

                        {/* Analysis Results (Molecules) */}
                        {msg.analysis && (
                            <div style={{ marginTop: '10px', display: 'flex', gap: '10px', flexWrap: 'wrap' }}>
                                {Array.isArray(msg.analysis) ? (
                                    msg.analysis.map((res: any, i: number) => (
                                        res.valid && res.structure_2d ? (
                                            <div key={i} style={{ width: '300px', height: '250px', position: 'relative' }}>
                                                <MoleculeViz structure={res.structure_2d} />
                                                <div style={{ position: 'absolute', top: 5, left: 5, background: 'rgba(0,0,0,0.7)', padding: '2px 6px', borderRadius: '4px', fontSize: '0.8em' }}>
                                                    {res.smiles}
                                                </div>
                                            </div>
                                        ) : null
                                    ))
                                ) : (
                                    msg.analysis.valid && msg.analysis.structure_2d && (
                                        <div style={{ width: '100%', maxWidth: '400px', height: '300px', position: 'relative' }}>
                                            <MoleculeViz structure={msg.analysis.structure_2d} width={400} height={300} />
                                            <div style={{ position: 'absolute', top: 5, left: 5, background: 'rgba(0,0,0,0.7)', padding: '2px 6px', borderRadius: '4px', fontSize: '0.8em' }}>
                                                {msg.analysis.smiles}
                                            </div>
                                        </div>
                                    )
                                )}
                            </div>
                        )}
                    </motion.div>
                ))}
                {isLoading && (
                    <motion.div
                        initial={{ opacity: 0 }}
                        animate={{ opacity: 1 }}
                        style={{ alignSelf: 'flex-start', color: 'var(--text-secondary)', fontStyle: 'italic' }}
                    >
                        Thinking...
                    </motion.div>
                )}
                <div ref={messagesEndRef} />
            </div>

            {/* Input Area */}
            <div className="input-area" style={{
                padding: '20px',
                background: 'rgba(15, 23, 42, 0.8)',
                borderTop: '1px solid rgba(255, 255, 255, 0.1)'
            }}>
                <div style={{ marginBottom: '10px', display: 'flex', gap: '10px' }}>
                    <button
                        className="btn-secondary"
                        style={{ fontSize: '0.8rem', padding: '4px 8px' }}
                        onClick={() => setInput("Compare the planarity and drug-likeness of these compounds: CCO, c1ccccc1, and a K3,3-like structure if possible (e.g. hypothetical cage).")}
                    >
                        Try: Compare Planarity
                    </button>
                    <button
                        className="btn-secondary"
                        style={{ fontSize: '0.8rem', padding: '4px 8px' }}
                        onClick={() => setInput("Analyze Aspirin (CC(=O)OC1=CC=CC=C1C(=O)O) and check its topology.")}
                    >
                        Try: Aspirin
                    </button>
                </div>
                <div style={{ display: 'flex', gap: '10px' }}>
                    <textarea
                        value={input}
                        onChange={(e) => setInput(e.target.value)}
                        onKeyDown={handleKeyDown}
                        placeholder="Type a message or paste SMILES (e.g., CCO)..."
                        style={{
                            flex: 1,
                            background: 'rgba(30, 41, 59, 0.5)',
                            border: '1px solid rgba(255, 255, 255, 0.1)',
                            borderRadius: '8px',
                            padding: '12px',
                            color: '#fff',
                            resize: 'none',
                            height: '50px',
                            outline: 'none'
                        }}
                    />
                    <button
                        onClick={handleSend}
                        disabled={isLoading || !input.trim()}
                        className="btn-primary"
                        style={{ height: '50px', padding: '0 24px' }}
                    >
                        Send
                    </button>
                </div>
            </div>
        </div>
    );
};

export default DrugDiscoveryChat;
