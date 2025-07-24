# 🚀 Prairie Genomics Suite R Shiny - Deployment Guide

## 🌟 **Quick Start for Sharing**

### Option 1: shinyapps.io (Easiest - FREE)

**1. Set up shinyapps.io account**
- Go to [shinyapps.io](https://www.shinyapps.io/)
- Sign up for free account (up to 5 apps, 25 hours/month)

**2. Install deployment tools**
```r
install.packages("rsconnect")
library(rsconnect)
```

**3. Connect your account**
- Go to your shinyapps.io dashboard
- Click "Tokens" → "Show" → copy token info
```r
rsconnect::setAccountInfo(
  name='YOUR_ACCOUNT_NAME',
  token='YOUR_TOKEN',
  secret='YOUR_SECRET'
)
```

**4. Deploy your app**
```r
# From your prairie-genomics-suite_shiny directory
rsconnect::deployApp(
  appDir = ".",
  appName = "prairie-genomics-suite",
  forceUpdate = TRUE
)
```

**5. Share the URL**: `https://your-account.shinyapps.io/prairie-genomics-suite/`

### Option 2: Local Network (Instant)

**For colleagues on your network:**
```r
# Start app with network access
shiny::runApp('app.R', host='0.0.0.0', port=3838)
```

**Find your IP address:**
- **Windows**: Open cmd → `ipconfig` → look for IPv4 Address
- **Mac**: System Preferences → Network → Advanced → TCP/IP
- **Linux**: Terminal → `hostname -I`

**Share URL**: `http://YOUR_IP_ADDRESS:3838`

Example: `http://192.168.1.100:3838`

### Option 3: GitHub + Binder (No deployment needed)

**1. Push to GitHub**
```bash
git init
git add .
git commit -m "Prairie Genomics Suite R Shiny app"
git remote add origin https://github.com/YOUR_USERNAME/prairie-genomics-suite-shiny.git
git push -u origin main
```

**2. Create runtime.txt file**
```
r-4.3-2023-10-01
```

**3. Share Binder link**: 
`https://mybinder.org/v2/gh/YOUR_USERNAME/prairie-genomics-suite-shiny/main?urlpath=shiny`

## 🐳 **Docker Deployment**

**For IT departments or advanced users:**

```bash
# Build container
docker build -t prairie-genomics-suite .

# Run locally
docker run -p 3838:3838 prairie-genomics-suite

# Access at http://localhost:3838/prairie-genomics-suite/
```

**Deploy to cloud:**
- **AWS**: ECS or EC2
- **Google Cloud**: Cloud Run
- **Azure**: Container Instances

## 💡 **Sharing Best Practices**

### For Research Groups
1. **Internal sharing**: Use Option 2 (local network)
2. **External collaborators**: Use Option 1 (shinyapps.io)
3. **Publication supplement**: Use Option 3 (GitHub + Binder)

### For IT Departments
1. **Enterprise deployment**: Use Docker + internal cloud
2. **High availability**: Load balancer + multiple containers
3. **Security**: VPN access + authentication

### For Conferences/Demos
1. **Live demo**: Local laptop + mobile hotspot
2. **Poster QR code**: Link to shinyapps.io deployment
3. **Workshop**: Pre-deployed on cloud with multiple instances

## 🔧 **Troubleshooting Deployment**

### shinyapps.io Issues
**Problem**: Deployment fails
```r
# Check logs
rsconnect::showLogs()

# Increase timeout
options(rsconnect.max.bundle.size = 5*1024^3)  # 5GB limit
```

**Problem**: App crashes on startup
- Check that all packages are in `install.R`
- Verify no hardcoded file paths
- Test locally first with `run_app.R`

### Local Network Issues
**Problem**: Others can't access
- Check firewall settings (allow port 3838)
- Verify IP address is correct
- Try different port: `port=8080`

**Problem**: App is slow
- Reduce sample data size
- Check available RAM
- Consider running on server instead of laptop

### Docker Issues
**Problem**: Build fails
```bash
# Check for errors
docker build -t prairie-genomics-suite . --no-cache

# Debug with interactive container
docker run -it rocker/shiny-verse:4.3.0 /bin/bash
```

## 📊 **Usage Analytics**

### shinyapps.io
- Built-in analytics dashboard
- Shows user sessions, geographic data
- Usage hours tracking

### Google Analytics (Optional)
Add to `app.R` head section:
```r
tags$head(
  HTML("<!-- Google Analytics code here -->")
)
```

## 🔒 **Security Considerations**

### Public Deployment
- ✅ No hardcoded passwords/tokens
- ✅ Input validation for file uploads
- ✅ Rate limiting for heavy computations
- ✅ Clear data handling policy

### Private Deployment
- 🔐 VPN or IP restriction
- 🔐 Authentication (Shiny Server Pro)
- 🔐 HTTPS certificates
- 🔐 Regular security updates

## 💰 **Cost Estimates**

### Free Options
- **shinyapps.io**: Free tier (5 apps, 25 hours/month)
- **Binder**: Free (but may be slow)
- **Local network**: Free (uses your hardware)

### Paid Options
- **shinyapps.io Pro**: $9/month (more apps, hours)
- **AWS EC2**: ~$10-50/month (depending on instance size)
- **Shiny Server Pro**: $995/year (enterprise features)

## 📞 **Support**

### Community Support
- **R Shiny Community**: [community.rstudio.com](https://community.rstudio.com/)
- **Stack Overflow**: Tag with `[r]` and `[shiny]`
- **GitHub Issues**: For app-specific problems

### Professional Support
- **RStudio**: Professional support plans
- **Consultants**: For custom deployment solutions

---

**🧬 Choose the deployment option that best fits your needs and technical comfort level!**