import uvicorn

if __name__ == "__main__":
    uvicorn.run("src.main:app", host="0.0.0.0", port=1010, log_level="info", reload=True, lifespan='on', timeout_keep_alive=60)
