/* Dark mode toggle with localStorage persistence */
$(document).ready(function () {
  var DARK_KEY = 'meDashboard_darkMode';

  // Restore preference on load
  if (localStorage.getItem(DARK_KEY) === '1') {
    $('body').addClass('dark-mode');
    updateToggleIcon(true);
  }

  $(document).on('click', '#dark-mode-toggle', function () {
    var isDark = $('body').toggleClass('dark-mode').hasClass('dark-mode');
    localStorage.setItem(DARK_KEY, isDark ? '1' : '0');
    updateToggleIcon(isDark);
  });

  function updateToggleIcon(isDark) {
    var btn = $('#dark-mode-toggle');
    if (isDark) {
      btn.html('<i class="fa fa-sun-o fa-fw"></i>');
      btn.attr('title', 'Switch to light mode');
    } else {
      btn.html('<i class="fa fa-moon-o fa-fw"></i>');
      btn.attr('title', 'Switch to dark mode');
    }
  }
});
