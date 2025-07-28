// Phase 1 - Modern JavaScript Interactions for Prairie Genomics Suite
// Enhanced UI interactions, real-time updates, and client-side optimization

// ===========================================
// INITIALIZATION AND GLOBAL UTILITIES
// ===========================================

// Initialize modern interactions when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
    initializeModernComponents();
    setupRealtimeUpdates();
    initializeAuthHandlers();
    initializeProgressTracking();
    console.log('🧬 Prairie Genomics Suite Phase 1 - Modern interactions initialized');
});

// Global utility functions
const PrairieGenomics = {
    // Enhanced notification system
    notify: function(message, type = 'info', duration = 5000) {
        const notification = document.createElement('div');
        notification.className = `alert-modern alert-modern-${type} notification-toast`;
        notification.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            z-index: 10000;
            min-width: 300px;
            animation: slideInRight 0.3s ease-out;
            box-shadow: var(--shadow-lg);
        `;
        
        notification.innerHTML = `
            <div class="alert-modern-icon">
                ${this.getNotificationIcon(type)}
            </div>
            <div class="alert-modern-content">
                <div class="alert-modern-message">${message}</div>
            </div>
            <button class="notification-close" onclick="this.parentElement.remove()" style="
                background: none; border: none; font-size: 1.25rem; cursor: pointer;
                color: inherit; opacity: 0.7; margin-left: auto;
            ">&times;</button>
        `;
        
        document.body.appendChild(notification);
        
        // Auto-remove after duration
        setTimeout(() => {
            if (notification.parentElement) {
                notification.style.animation = 'slideOutRight 0.3s ease-in';
                setTimeout(() => notification.remove(), 300);
            }
        }, duration);
    },
    
    getNotificationIcon: function(type) {
        const icons = {
            success: '✅',
            error: '❌',
            warning: '⚠️',
            info: 'ℹ️'
        };
        return icons[type] || icons.info;
    },
    
    // Enhanced loading states
    showLoading: function(element, message = 'Loading...') {
        if (typeof element === 'string') {
            element = document.querySelector(element);
        }
        
        if (!element) return;
        
        element.classList.add('loading-state');
        element.innerHTML = `
            <div class="loading-modern">
                <div class="loading-modern-spinner"></div>
                <div class="loading-modern-text">${message}</div>
            </div>
        `;
    },
    
    hideLoading: function(element, originalContent = '') {
        if (typeof element === 'string') {
            element = document.querySelector(element);
        }
        
        if (!element) return;
        
        element.classList.remove('loading-state');
        element.innerHTML = originalContent;
    },
    
    // Progress tracking utilities
    updateProgress: function(elementId, percentage, message = '') {
        const progressContainer = document.getElementById(elementId);
        if (!progressContainer) return;
        
        const progressBar = progressContainer.querySelector('.progress-modern-bar');
        const progressLabel = progressContainer.querySelector('.progress-modern-label span:last-child');
        const progressMessage = progressContainer.querySelector('.progress-modern-message');
        
        if (progressBar) {
            progressBar.style.width = `${Math.min(100, Math.max(0, percentage))}%`;
        }
        
        if (progressLabel) {
            progressLabel.textContent = `${Math.round(percentage)}%`;
        }
        
        if (progressMessage && message) {
            progressMessage.textContent = message;
        }
    }
};

// ===========================================
// MODERN COMPONENT INITIALIZATION
// ===========================================

function initializeModernComponents() {
    // Initialize tooltips
    initializeTooltips();
    
    // Initialize smooth scrolling
    initializeSmoothScrolling();
    
    // Initialize enhanced form interactions
    initializeFormEnhancements();
    
    // Initialize table enhancements
    initializeTableEnhancements();
    
    // Initialize card animations
    initializeCardAnimations();
}

function initializeTooltips() {
    const tooltipElements = document.querySelectorAll('[data-tooltip]');
    
    tooltipElements.forEach(element => {
        element.classList.add('tooltip-modern');
        
        // Add hover delay for better UX
        let hoverTimer;
        
        element.addEventListener('mouseenter', function() {
            hoverTimer = setTimeout(() => {
                this.classList.add('tooltip-visible');
            }, 500);
        });
        
        element.addEventListener('mouseleave', function() {
            clearTimeout(hoverTimer);
            this.classList.remove('tooltip-visible');
        });
    });
}

function initializeSmoothScrolling() {
    // Add smooth scrolling to anchor links
    const anchorLinks = document.querySelectorAll('a[href^="#"]');
    
    anchorLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            const href = this.getAttribute('href');
            const target = document.querySelector(href);
            
            if (target) {
                e.preventDefault();
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });
}

function initializeFormEnhancements() {
    // Add floating labels and validation styling
    const formInputs = document.querySelectorAll('.form-input-modern, .form-select-modern');
    
    formInputs.forEach(input => {
        // Add focus animations
        input.addEventListener('focus', function() {
            this.parentElement.classList.add('form-group-focused');
        });
        
        input.addEventListener('blur', function() {
            this.parentElement.classList.remove('form-group-focused');
            
            // Add validation styling
            if (this.validity && !this.validity.valid) {
                this.classList.add('form-input-error');
            } else {
                this.classList.remove('form-input-error');
            }
        });
        
        // Real-time validation for enhanced UX
        input.addEventListener('input', function() {
            if (this.classList.contains('form-input-error') && this.validity.valid) {
                this.classList.remove('form-input-error');
            }
        });
    });
}

function initializeTableEnhancements() {
    const tables = document.querySelectorAll('.table-modern');
    
    tables.forEach(table => {
        // Add sortable headers
        const headers = table.querySelectorAll('th[data-sortable]');
        
        headers.forEach(header => {
            header.style.cursor = 'pointer';
            header.innerHTML += ' <span class="sort-indicator">↕️</span>';
            
            header.addEventListener('click', function() {
                sortTable(table, this);
            });
        });
        
        // Add row selection
        const rows = table.querySelectorAll('tbody tr');
        rows.forEach(row => {
            row.addEventListener('click', function() {
                // Toggle row selection
                this.classList.toggle('row-selected');
                
                // Dispatch custom event for row selection
                const event = new CustomEvent('rowSelected', {
                    detail: {
                        row: this,
                        data: this.dataset
                    }
                });
                table.dispatchEvent(event);
            });
        });
    });
}

function initializeCardAnimations() {
    const cards = document.querySelectorAll('.card-modern');
    
    // Intersection Observer for scroll animations
    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };
    
    const observer = new IntersectionObserver((entries) => {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('card-animated');
            }
        });
    }, observerOptions);
    
    cards.forEach(card => {
        observer.observe(card);
    });
}

// ===========================================
// REAL-TIME UPDATES AND WEBSOCKETS
// ===========================================

function setupRealtimeUpdates() {
    // Set up Shiny custom message handlers
    if (window.Shiny) {
        // Progress update handler
        Shiny.addCustomMessageHandler('updateProgress', function(data) {
            PrairieGenomics.updateProgress(data.elementId, data.percentage, data.message);
        });
        
        // Notification handler
        Shiny.addCustomMessageHandler('showNotification', function(data) {
            PrairieGenomics.notify(data.message, data.type, data.duration);
        });
        
        // Loading state handler
        Shiny.addCustomMessageHandler('setLoadingState', function(data) {
            if (data.show) {
                PrairieGenomics.showLoading(data.elementId, data.message);
            } else {
                PrairieGenomics.hideLoading(data.elementId, data.content);
            }
        });
        
        // Analysis progress handler
        Shiny.addCustomMessageHandler('analysisProgress', function(data) {
            updateAnalysisProgress(data);
        });
        
        // Real-time results update
        Shiny.addCustomMessageHandler('updateResults', function(data) {
            updateResultsDisplay(data);
        });
    }
}

function updateAnalysisProgress(data) {
    const progressContainer = document.getElementById('analysis-progress-container');
    if (!progressContainer) return;
    
    // Update main progress bar
    PrairieGenomics.updateProgress('main-progress', data.percentage, data.message);
    
    // Update step indicators
    const steps = document.querySelectorAll('.step-indicator');
    steps.forEach((step, index) => {
        if (index < data.currentStep) {
            step.classList.add('step-completed');
            step.classList.remove('step-active', 'step-pending');
        } else if (index === data.currentStep) {
            step.classList.add('step-active');
            step.classList.remove('step-completed', 'step-pending');
        } else {
            step.classList.add('step-pending');
            step.classList.remove('step-completed', 'step-active');
        }
    });
    
    // Update detailed status
    const statusElement = document.getElementById('analysis-status-detail');
    if (statusElement) {
        statusElement.innerHTML = `
            <div class="status-item">
                <strong>Current Step:</strong> ${data.stepName}
            </div>
            <div class="status-item">
                <strong>Time Elapsed:</strong> ${formatDuration(data.timeElapsed)}
            </div>
            <div class="status-item">
                <strong>ETA:</strong> ${data.eta || 'Calculating...'}
            </div>
        `;
    }
}

function updateResultsDisplay(data) {
    const resultsContainer = document.getElementById('results-container');
    if (!resultsContainer) return;
    
    // Smooth transition for results update
    resultsContainer.style.opacity = '0.5';
    
    setTimeout(() => {
        // Update results content
        if (data.html) {
            resultsContainer.innerHTML = data.html;
        }
        
        // Update summary statistics
        if (data.stats) {
            updateSummaryStats(data.stats);
        }
        
        // Reinitialize components for new content
        initializeModernComponents();
        
        // Fade back in
        resultsContainer.style.opacity = '1';
        
        // Show success notification
        PrairieGenomics.notify('Results updated successfully', 'success');
    }, 300);
}

// ===========================================
// FIREBASE AUTHENTICATION HANDLERS
// ===========================================

function initializeAuthHandlers() {
    // Set up authentication form handlers
    const signInForm = document.getElementById('signin-form');
    const signUpForm = document.getElementById('signup-form');
    
    if (signInForm) {
        signInForm.addEventListener('submit', handleSignIn);
    }
    
    if (signUpForm) {
        signUpForm.addEventListener('submit', handleSignUp);
    }
    
    // Set up auth state listeners
    setupAuthStateListeners();
    
    // Set up custom message handlers for auth
    if (window.Shiny) {
        Shiny.addCustomMessageHandler('store_auth_token', function(data) {
            localStorage.setItem('prairie_auth_token', data.token);
            localStorage.setItem('prairie_user_data', JSON.stringify(data.user));
            updateAuthUI(data.user);
        });
        
        Shiny.addCustomMessageHandler('clear_auth_token', function(data) {
            localStorage.removeItem('prairie_auth_token');
            localStorage.removeItem('prairie_user_data');
            updateAuthUI(null);
        });
    }
}

function handleSignIn(event) {
    event.preventDefault();
    
    const formData = new FormData(event.target);
    const email = formData.get('email');
    const password = formData.get('password');
    
    // Show loading state
    const submitButton = event.target.querySelector('button[type="submit"]');
    const originalText = submitButton.textContent;
    submitButton.classList.add('btn-modern-loading');
    submitButton.disabled = true;
    
    // Send to Shiny for processing
    if (window.Shiny) {
        Shiny.setInputValue('auth_signin', {
            email: email,
            password: password,
            timestamp: Date.now()
        });
    }
    
    // Reset button after timeout
    setTimeout(() => {
        submitButton.classList.remove('btn-modern-loading');
        submitButton.disabled = false;
        submitButton.textContent = originalText;
    }, 5000);
}

function handleSignUp(event) {
    event.preventDefault();
    
    const formData = new FormData(event.target);
    const email = formData.get('email');
    const password = formData.get('password');
    const confirmPassword = formData.get('confirm_password');
    const displayName = formData.get('display_name');
    
    // Client-side validation
    if (password !== confirmPassword) {
        PrairieGenomics.notify('Passwords do not match', 'error');
        return;
    }
    
    if (password.length < 6) {
        PrairieGenomics.notify('Password must be at least 6 characters', 'error');
        return;
    }
    
    // Show loading state
    const submitButton = event.target.querySelector('button[type="submit"]');
    submitButton.classList.add('btn-modern-loading');
    submitButton.disabled = true;
    
    // Send to Shiny for processing
    if (window.Shiny) {
        Shiny.setInputValue('auth_signup', {
            email: email,
            password: password,
            displayName: displayName,
            timestamp: Date.now()
        });
    }
}

function setupAuthStateListeners() {
    // Check for existing auth token on page load
    const token = localStorage.getItem('prairie_auth_token');
    const userData = localStorage.getItem('prairie_user_data');
    
    if (token && userData) {
        try {
            const user = JSON.parse(userData);
            updateAuthUI(user);
            
            // Verify token with server
            if (window.Shiny) {
                Shiny.setInputValue('verify_auth_token', {
                    token: token,
                    timestamp: Date.now()
                });
            }
        } catch (e) {
            // Clear invalid data
            localStorage.removeItem('prairie_auth_token');
            localStorage.removeItem('prairie_user_data');
        }
    }
}

function updateAuthUI(user) {
    const authSection = document.getElementById('auth-section');
    const userSection = document.getElementById('user-section');
    const userDisplay = document.getElementById('user-display');
    
    if (user) {
        // User is authenticated
        if (authSection) authSection.style.display = 'none';
        if (userSection) userSection.style.display = 'block';
        if (userDisplay) {
            userDisplay.innerHTML = `
                <div class="user-avatar">
                    ${user.displayName.charAt(0).toUpperCase()}
                </div>
                <div class="user-info">
                    <div class="user-name">${user.displayName}</div>
                    <div class="user-email">${user.email}</div>
                </div>
            `;
        }
        
        // Enable authenticated features
        enableAuthenticatedFeatures();
        
    } else {
        // User is not authenticated
        if (authSection) authSection.style.display = 'block';
        if (userSection) userSection.style.display = 'none';
        
        // Disable authenticated features
        disableAuthenticatedFeatures();
    }
}

// ===========================================
// PROGRESS TRACKING AND ASYNC OPERATIONS
// ===========================================

function initializeProgressTracking() {
    // Set up progress tracking for async operations
    if (window.Shiny) {
        // DESeq2 analysis progress
        Shiny.addCustomMessageHandler('deseq2_progress', function(data) {
            updateAnalysisProgress({
                percentage: data.percentage,
                message: data.message,
                currentStep: data.step,
                stepName: data.stepName,
                timeElapsed: data.timeElapsed,
                eta: data.eta
            });
        });
        
        // Gene conversion progress
        Shiny.addCustomMessageHandler('gene_conversion_progress', function(data) {
            PrairieGenomics.updateProgress('gene-conversion-progress', data.percentage, data.message);
        });
        
        // Pathway analysis progress
        Shiny.addCustomMessageHandler('pathway_progress', function(data) {
            PrairieGenomics.updateProgress('pathway-analysis-progress', data.percentage, data.message);
        });
    }
}

// ===========================================
// UTILITY FUNCTIONS
// ===========================================

function formatDuration(milliseconds) {
    const seconds = Math.floor(milliseconds / 1000);
    const minutes = Math.floor(seconds / 60);
    const hours = Math.floor(minutes / 60);
    
    if (hours > 0) {
        return `${hours}h ${minutes % 60}m ${seconds % 60}s`;
    } else if (minutes > 0) {
        return `${minutes}m ${seconds % 60}s`;
    } else {
        return `${seconds}s`;
    }
}

function sortTable(table, header) {
    const headerIndex = Array.from(header.parentElement.children).indexOf(header);
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    
    // Determine sort direction
    const currentDirection = header.dataset.sortDirection || 'asc';
    const newDirection = currentDirection === 'asc' ? 'desc' : 'asc';
    header.dataset.sortDirection = newDirection;
    
    // Update sort indicators
    table.querySelectorAll('th .sort-indicator').forEach(indicator => {
        indicator.textContent = '↕️';
    });
    header.querySelector('.sort-indicator').textContent = newDirection === 'asc' ? '↑' : '↓';
    
    // Sort rows
    rows.sort((a, b) => {
        const aValue = a.children[headerIndex].textContent.trim();
        const bValue = b.children[headerIndex].textContent.trim();
        
        // Try numeric comparison first
        const aNum = parseFloat(aValue);
        const bNum = parseFloat(bValue);
        
        if (!isNaN(aNum) && !isNaN(bNum)) {
            return newDirection === 'asc' ? aNum - bNum : bNum - aNum;
        } else {
            // String comparison
            return newDirection === 'asc' 
                ? aValue.localeCompare(bValue)
                : bValue.localeCompare(aValue);
        }
    });
    
    // Rebuild tbody
    rows.forEach(row => tbody.appendChild(row));
}

function updateSummaryStats(stats) {
    const statsContainer = document.getElementById('summary-stats');
    if (!statsContainer) return;
    
    const statsHtml = Object.entries(stats).map(([key, value]) => `
        <div class="stat-item">
            <div class="stat-value">${value.toLocaleString()}</div>
            <div class="stat-label">${key.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase())}</div>
        </div>
    `).join('');
    
    statsContainer.innerHTML = statsHtml;
}

function enableAuthenticatedFeatures() {
    const authRequiredElements = document.querySelectorAll('[data-auth-required]');
    authRequiredElements.forEach(element => {
        element.style.display = '';
        element.disabled = false;
    });
}

function disableAuthenticatedFeatures() {
    const authRequiredElements = document.querySelectorAll('[data-auth-required]');
    authRequiredElements.forEach(element => {
        element.style.display = 'none';
        element.disabled = true;
    });
}

// ===========================================
// CSS ANIMATIONS
// ===========================================

// Add CSS for animations
const styleSheet = document.createElement('style');
styleSheet.textContent = `
    @keyframes slideInRight {
        from {
            transform: translateX(100%);
            opacity: 0;
        }
        to {
            transform: translateX(0);
            opacity: 1;
        }
    }
    
    @keyframes slideOutRight {
        from {
            transform: translateX(0);
            opacity: 1;
        }
        to {
            transform: translateX(100%);
            opacity: 0;
        }
    }
    
    .card-animated {
        animation: fadeInUp 0.6s ease-out;
    }
    
    @keyframes fadeInUp {
        from {
            transform: translateY(30px);
            opacity: 0;
        }
        to {
            transform: translateY(0);
            opacity: 1;
        }
    }
    
    .form-group-focused .form-label-modern {
        color: var(--primary-600);
        transform: translateY(-2px);
        transition: all var(--transition-fast);
    }
    
    .form-input-error {
        border-color: var(--error-500) !important;
        box-shadow: 0 0 0 3px rgba(239, 68, 68, 0.1) !important;
    }
    
    .row-selected {
        background: var(--primary-50) !important;
        border-left: 3px solid var(--primary-500);
    }
    
    .step-indicator {
        display: flex;
        align-items: center;
        padding: var(--space-3);
        border-radius: var(--rounded-lg);
        margin-bottom: var(--space-2);
        transition: all var(--transition-fast);
    }
    
    .step-completed {
        background: var(--success-50);
        color: var(--success-700);
    }
    
    .step-active {
        background: var(--primary-50);
        color: var(--primary-700);
        border-left: 3px solid var(--primary-500);
    }
    
    .step-pending {
        background: var(--gray-50);
        color: var(--gray-500);
    }
    
    .user-avatar {
        width: 2.5rem;
        height: 2.5rem;
        border-radius: var(--rounded-full);
        background: var(--primary-500);
        color: white;
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: 600;
        margin-right: var(--space-3);
    }
    
    .stat-item {
        text-align: center;
        padding: var(--space-4);
        border-radius: var(--rounded-lg);
        background: var(--gray-50);
    }
    
    .stat-value {
        font-size: 2rem;
        font-weight: 700;
        color: var(--primary-600);
        margin-bottom: var(--space-1);
    }
    
    .stat-label {
        font-size: 0.875rem;
        color: var(--gray-600);
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
`;

document.head.appendChild(styleSheet);

// Export for global access
window.PrairieGenomics = PrairieGenomics;