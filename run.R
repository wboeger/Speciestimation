# run.R -- container entrypoint. Reads Railway's $PORT (falls back to 8080
# for local `docker run`) and serves the app directory on 0.0.0.0.
port <- as.integer(Sys.getenv("PORT", "8080"))
shiny::runApp(appDir = ".", host = "0.0.0.0", port = port)
