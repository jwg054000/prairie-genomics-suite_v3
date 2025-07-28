/**
 * Prairie Genomics Suite - Phase 2 Real-time Client
 * Handles real-time updates, notifications, and WebSocket-style communication
 */

// Global namespace for Prairie Genomics real-time features
window.PrairieRealtimeClient = (function() {
  'use strict';

  // Configuration
  const config = {
    notificationDuration: 5000,
    maxNotifications: 5,
    updateInterval: 1000,
    animationDuration: 300
  };

  // State management
  const state = {
    activeNotifications: new Map(),
    activeProcesses: new Map(),
    systemMetrics: {},
    isInitialized: false
  };

  // ===========================================
  // INITIALIZATION
  // ===========================================

  function initialize() {
    if (state.isInitialized) return;

    console.log('🚀 Initializing Prairie Genomics Real-time Client');

    // Create notification container if it doesn't exist
    createNotificationContainer();

    // Set up Shiny message handlers
    setupShinyMessageHandlers();

    // Initialize UI components
    initializeRealtimeUI();

    // Start heartbeat
    startHeartbeat();

    state.isInitialized = true;
    console.log('✅ Real-time client initialized successfully');
  }

  // ===========================================
  // NOTIFICATION SYSTEM
  // ===========================================

  function createNotificationContainer() {
    let container = document.getElementById('realtime-notifications');
    if (!container) {
      container = document.createElement('div');
      container.id = 'realtime-notifications';
      container.className = 'realtime-notifications';
      container.style.cssText = `
        position: fixed;
        top: 20px;
        right: 20px;
        z-index: 9999;
        width: 350px;
        pointer-events: none;
      `;
      document.body.appendChild(container);
    }
  }

  function showNotification(data) {
    const { id, type, title, message, duration = config.notificationDuration, action } = data;

    // Remove old notification if exists
    if (state.activeNotifications.has(id)) {
      removeNotification(id);
    }

    // Create notification element
    const notification = createNotificationElement(id, type, title, message, action);
    
    // Add to DOM
    const container = document.getElementById('realtime-notifications');
    container.appendChild(notification);

    // Add to state
    state.activeNotifications.set(id, {
      element: notification,
      timer: setTimeout(() => removeNotification(id), duration)
    });

    // Animate in
    setTimeout(() => {
      notification.classList.add('show');
    }, 10);

    // Limit number of notifications
    while (state.activeNotifications.size > config.maxNotifications) {
      const oldestId = state.activeNotifications.keys().next().value;
      removeNotification(oldestId);
    }

    console.log(`📢 Notification shown: ${type} - ${title}`);
  }

  function createNotificationElement(id, type, title, message, action) {
    const notification = document.createElement('div');
    notification.id = `notification-${id}`;
    notification.className = `realtime-notification notification-${type}`;
    
    // Get type-specific styling and icon
    const typeConfig = getNotificationTypeConfig(type);
    
    notification.style.cssText = `
      background: ${typeConfig.background};
      border-left: 4px solid ${typeConfig.borderColor};
      color: ${typeConfig.textColor};
      margin-bottom: 12px;
      padding: 16px;
      border-radius: 8px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.15);
      transform: translateX(100%);
      transition: all ${config.animationDuration}ms cubic-bezier(0.4, 0, 0.2, 1);
      opacity: 0;
      pointer-events: auto;
      position: relative;
      overflow: hidden;
    `;

    // Progress bar for duration
    const progressBar = document.createElement('div');
    progressBar.className = 'notification-progress';
    progressBar.style.cssText = `
      position: absolute;
      bottom: 0;
      left: 0;
      height: 3px;
      background: ${typeConfig.progressColor};
      width: 100%;
      transform-origin: left;
      animation: notificationProgress ${config.notificationDuration}ms linear forwards;
    `;

    notification.innerHTML = `
      <div style="display: flex; align-items: flex-start;">
        <div style="margin-right: 12px; font-size: 20px;">
          ${typeConfig.icon}
        </div>
        <div style="flex: 1;">
          <div style="font-weight: 600; margin-bottom: 4px; font-size: 15px;">
            ${escapeHtml(title)}
          </div>
          <div style="font-size: 14px; opacity: 0.9; line-height: 1.4;">
            ${escapeHtml(message)}
          </div>
          ${action ? `
            <div style="margin-top: 12px;">
              <button class="notification-action-btn" 
                      style="background: ${typeConfig.buttonColor}; color: white; 
                             border: none; padding: 6px 12px; border-radius: 4px; 
                             font-size: 13px; cursor: pointer;">
                ${escapeHtml(action.label || 'Action')}
              </button>
            </div>
          ` : ''}
        </div>
        <button class="notification-close" 
                style="background: none; border: none; color: inherit; 
                       opacity: 0.7; cursor: pointer; font-size: 18px; 
                       padding: 0; margin-left: 8px;"
                onclick="PrairieRealtimeClient.removeNotification('${id}')">
          ×
        </button>
      </div>
    `;

    notification.appendChild(progressBar);

    // Add show class for animation
    notification.classList.add('notification-enter');

    return notification;
  }

  function getNotificationTypeConfig(type) {
    const configs = {
      success: {
        icon: '✅',
        background: '#f0f9ff',
        borderColor: '#10b981',
        textColor: '#064e3b',
        progressColor: '#10b981',
        buttonColor: '#10b981'
      },
      error: {
        icon: '❌',
        background: '#fef2f2',
        borderColor: '#ef4444',
        textColor: '#7f1d1d',
        progressColor: '#ef4444',
        buttonColor: '#ef4444'
      },
      warning: {
        icon: '⚠️',
        background: '#fffbeb',
        borderColor: '#f59e0b',
        textColor: '#78350f',
        progressColor: '#f59e0b',
        buttonColor: '#f59e0b'
      },
      info: {
        icon: 'ℹ️',
        background: '#f0f9ff',
        borderColor: '#3b82f6',
        textColor: '#1e3a8a',
        progressColor: '#3b82f6',
        buttonColor: '#3b82f6'
      }
    };

    return configs[type] || configs.info;
  }

  function removeNotification(id) {
    const notificationData = state.activeNotifications.get(id);
    if (!notificationData) return;

    const { element, timer } = notificationData;

    // Clear timer
    if (timer) clearTimeout(timer);

    // Animate out
    element.style.transform = 'translateX(100%)';
    element.style.opacity = '0';

    // Remove from DOM after animation
    setTimeout(() => {
      if (element.parentNode) {
        element.parentNode.removeChild(element);
      }
    }, config.animationDuration);

    // Remove from state
    state.activeNotifications.delete(id);

    console.log(`🗑️ Notification removed: ${id}`);
  }

  function clearAllNotifications() {
    state.activeNotifications.forEach((_, id) => {
      removeNotification(id);
    });
    console.log('🧹 All notifications cleared');
  }

  // ===========================================
  // PROGRESS TRACKING
  // ===========================================

  function updateProcessProgress(data) {
    const { process_id, progress, message, elapsed_time, eta, current_step, total_steps } = data;

    // Update or create progress element
    let progressElement = document.getElementById(`progress-${process_id}`);
    if (!progressElement) {
      progressElement = createProgressElement(process_id);
    }

    // Update progress bar
    const progressBar = progressElement.querySelector('.process-progress-bar');
    const progressText = progressElement.querySelector('.process-progress-text');
    const progressMessage = progressElement.querySelector('.process-progress-message');
    const progressDetails = progressElement.querySelector('.process-progress-details');

    if (progressBar) {
      progressBar.style.width = `${progress}%`;
      progressBar.setAttribute('aria-valuenow', progress);
    }

    if (progressText) {
      progressText.textContent = `${Math.round(progress)}%`;
    }

    if (progressMessage) {
      progressMessage.textContent = message || 'Processing...';
    }

    if (progressDetails) {
      let detailsText = '';
      if (elapsed_time !== undefined) {
        detailsText += `Elapsed: ${formatTime(elapsed_time)}`;
      }
      if (eta !== undefined && eta > 0) {
        detailsText += ` | ETA: ${formatTime(eta)}`;
      }
      if (current_step && total_steps) {
        detailsText += ` | Step ${current_step}/${total_steps}`;
      }
      progressDetails.textContent = detailsText;
    }

    // Update state
    state.activeProcesses.set(process_id, data);

    console.log(`📊 Progress updated: ${process_id} - ${progress}%`);
  }

  function createProgressElement(process_id) {
    const container = document.getElementById('active-processes-list') || document.body;
    
    const progressElement = document.createElement('div');
    progressElement.id = `progress-${process_id}`;
    progressElement.className = 'process-progress-container';
    progressElement.style.cssText = `
      background: white;
      border: 1px solid #e5e7eb;
      border-radius: 8px;
      padding: 16px;
      margin-bottom: 12px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    `;

    progressElement.innerHTML = `
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
        <div class="process-progress-message" style="font-weight: 500; color: #374151;">
          Processing...
        </div>
        <div class="process-progress-text" style="font-weight: 600; color: #3b82f6;">
          0%
        </div>
      </div>
      <div class="progress-modern" style="background: #f3f4f6; height: 8px; border-radius: 4px; overflow: hidden;">
        <div class="process-progress-bar" 
             style="background: linear-gradient(90deg, #3b82f6, #8b5cf6); 
                    height: 100%; width: 0%; transition: width 0.3s ease;
                    border-radius: 4px;"
             role="progressbar" aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">
        </div>
      </div>
      <div class="process-progress-details" 
           style="font-size: 12px; color: #6b7280; margin-top: 6px;">
      </div>
    `;

    container.appendChild(progressElement);
    return progressElement;
  }

  function completeProcess(data) {
    const { process_id, total_duration, final_message } = data;

    const progressElement = document.getElementById(`progress-${process_id}`);
    if (progressElement) {
      // Update to completion state
      const progressBar = progressElement.querySelector('.process-progress-bar');
      const progressText = progressElement.querySelector('.process-progress-text');
      const progressMessage = progressElement.querySelector('.process-progress-message');
      const progressDetails = progressElement.querySelector('.process-progress-details');

      if (progressBar) {
        progressBar.style.width = '100%';
        progressBar.style.background = 'linear-gradient(90deg, #10b981, #059669)';
      }

      if (progressText) {
        progressText.textContent = '✓ Complete';
        progressText.style.color = '#10b981';
      }

      if (progressMessage) {
        progressMessage.textContent = final_message || 'Complete';
        progressMessage.style.color = '#10b981';
      }

      if (progressDetails) {
        progressDetails.textContent = `Completed in ${formatTime(total_duration)}`;
      }

      // Remove after delay
      setTimeout(() => {
        progressElement.style.opacity = '0';
        progressElement.style.transform = 'translateY(-10px)';
        setTimeout(() => {
          if (progressElement.parentNode) {
            progressElement.parentNode.removeChild(progressElement);
          }
        }, config.animationDuration);
      }, 3000);
    }

    // Remove from active processes
    state.activeProcesses.delete(process_id);

    console.log(`✅ Process completed: ${process_id} (${total_duration}s)`);
  }

  function failProcess(data) {
    const { process_id, error_message, total_duration } = data;

    const progressElement = document.getElementById(`progress-${process_id}`);
    if (progressElement) {
      // Update to error state
      const progressBar = progressElement.querySelector('.process-progress-bar');
      const progressText = progressElement.querySelector('.process-progress-text');
      const progressMessage = progressElement.querySelector('.process-progress-message');
      const progressDetails = progressElement.querySelector('.process-progress-details');

      if (progressBar) {
        progressBar.style.background = 'linear-gradient(90deg, #ef4444, #dc2626)';
      }

      if (progressText) {
        progressText.textContent = '✗ Failed';
        progressText.style.color = '#ef4444';
      }

      if (progressMessage) {
        progressMessage.textContent = error_message || 'Process failed';
        progressMessage.style.color = '#ef4444';
      }

      if (progressDetails) {
        progressDetails.textContent = `Failed after ${formatTime(total_duration)}`;
      }

      // Keep visible longer for errors
      setTimeout(() => {
        progressElement.style.opacity = '0.7';
        progressElement.style.border = '1px solid #fee2e2';
      }, 8000);
    }

    // Remove from active processes
    state.activeProcesses.delete(process_id);

    console.log(`❌ Process failed: ${process_id} - ${error_message}`);
  }

  // ===========================================
  // SYSTEM METRICS
  // ===========================================

  function updateSystemMetrics(data) {
    const { memory_usage, memory_percent, server_load, active_sessions, last_update } = data;

    // Update state
    state.systemMetrics = data;

    // Update UI elements
    updateElement('memory-usage', `${memory_percent}%`);
    updateElement('last-update', last_update);
    updateElement('active-count', state.activeProcesses.size.toString());

    // Update status indicator
    const statusIndicator = document.getElementById('global-status');
    if (statusIndicator) {
      if (state.activeProcesses.size > 0) {
        statusIndicator.textContent = '●';
        statusIndicator.style.color = '#f59e0b'; // Orange for active
      } else if (memory_percent > 80) {
        statusIndicator.textContent = '●';
        statusIndicator.style.color = '#ef4444'; // Red for high memory
      } else {
        statusIndicator.textContent = '●';
        statusIndicator.style.color = '#10b981'; // Green for idle
      }
    }

    console.log(`📊 System metrics updated: Memory ${memory_percent}%, Load ${server_load}%`);
  }

  // ===========================================
  // SHINY MESSAGE HANDLERS
  // ===========================================

  function setupShinyMessageHandlers() {
    if (typeof Shiny === 'undefined') {
      console.warn('⚠️ Shiny not found, message handlers not registered');
      return;
    }

    // Notification handlers
    Shiny.addCustomMessageHandler('showRealtimeNotification', showNotification);
    Shiny.addCustomMessageHandler('clearAllNotifications', clearAllNotifications);

    // Progress handlers
    Shiny.addCustomMessageHandler('updateProcessProgress', updateProcessProgress);
    Shiny.addCustomMessageHandler('completeProcess', completeProcess);
    Shiny.addCustomMessageHandler('failProcess', failProcess);
    Shiny.addCustomMessageHandler('clearProcessHistory', clearProcessHistory);

    // System handlers
    Shiny.addCustomMessageHandler('updateSystemMetrics', updateSystemMetrics);
    Shiny.addCustomMessageHandler('updateGlobalProgress', updateGlobalProgress);

    console.log('📡 Shiny message handlers registered');
  }

  // ===========================================
  // UTILITY FUNCTIONS
  // ===========================================

  function updateElement(id, content) {
    const element = document.getElementById(id);
    if (element) {
      element.textContent = content;
    }
  }

  function formatTime(seconds) {
    if (seconds < 60) {
      return `${Math.round(seconds)}s`;
    } else if (seconds < 3600) {
      const minutes = Math.floor(seconds / 60);
      const remainingSeconds = Math.round(seconds % 60);
      return `${minutes}m ${remainingSeconds}s`;
    } else {
      const hours = Math.floor(seconds / 3600);
      const remainingMinutes = Math.floor((seconds % 3600) / 60);
      return `${hours}h ${remainingMinutes}m`;
    }
  }

  function escapeHtml(text) {
    const div = document.createElement('div');
    div.textContent = text;
    return div.innerHTML;
  }

  function clearProcessHistory() {
    // Remove all completed/failed process elements
    const processElements = document.querySelectorAll('[id^="progress-"]');
    processElements.forEach(element => {
      if (element.parentNode) {
        element.parentNode.removeChild(element);
      }
    });
    console.log('🧹 Process history cleared');
  }

  function updateGlobalProgress(data) {
    const { status, active_count, total_count } = data;
    
    // Show/hide active processes panel
    const panel = document.getElementById('active-processes-panel');
    if (panel) {
      panel.style.display = active_count > 0 ? 'block' : 'none';
    }
  }

  function initializeRealtimeUI() {
    // Add CSS animations
    const style = document.createElement('style');
    style.textContent = `
      @keyframes notificationProgress {
        from { transform: scaleX(1); }
        to { transform: scaleX(0); }
      }

      .realtime-notification.show {
        transform: translateX(0) !important;
        opacity: 1 !important;
      }

      .realtime-notification:hover .notification-progress {
        animation-play-state: paused;
      }

      .notification-action-btn:hover {
        opacity: 0.9;
        transform: translateY(-1px);
      }
    `;
    document.head.appendChild(style);
  }

  function startHeartbeat() {
    // Send periodic heartbeat to maintain connection
    setInterval(() => {
      if (typeof Shiny !== 'undefined' && Shiny.shinyapp) {
        // Send heartbeat message to server
        try {
          Shiny.setInputValue('realtime_heartbeat', new Date().getTime(), {priority: 'event'});
        } catch (e) {
          console.warn('⚠️ Heartbeat failed:', e.message);
        }
      }
    }, 30000); // Every 30 seconds
  }

  // ===========================================
  // PUBLIC API
  // ===========================================

  return {
    initialize: initialize,
    showNotification: showNotification,
    removeNotification: removeNotification,
    clearAllNotifications: clearAllNotifications,
    updateProcessProgress: updateProcessProgress,
    completeProcess: completeProcess,
    failProcess: failProcess,
    updateSystemMetrics: updateSystemMetrics,
    getState: () => ({ ...state })
  };

})();

// Auto-initialize when DOM is ready
document.addEventListener('DOMContentLoaded', function() {
  PrairieRealtimeClient.initialize();
});

// Initialize immediately if DOM already loaded
if (document.readyState === 'loading') {
  // DOM not ready, event listener will handle it
} else {
  // DOM is ready
  PrairieRealtimeClient.initialize();
}