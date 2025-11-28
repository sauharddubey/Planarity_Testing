#!/bin/bash

# Start Backend
cd backend
uvicorn main:app --reload --port 8000 &
BACKEND_PID=$!

# Start Frontend
cd ../frontend
npm run dev -- --port 3000 &
FRONTEND_PID=$!

echo "Backend running on port 8000 (PID: $BACKEND_PID)"
echo "Frontend running on port 3000 (PID: $FRONTEND_PID)"

trap "kill $BACKEND_PID $FRONTEND_PID" EXIT

wait
